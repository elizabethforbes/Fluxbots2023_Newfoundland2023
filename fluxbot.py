import pandas as pd
import numpy as np
import plotly.graph_objs as go
# import chart_studio.plotly as py
import plotly
# from sklearn.linear_model import LinearRegression
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.pipeline import make_pipeline
from scipy import optimize
import math

treatments = {
    'MWC': 'Full Exclusion',
    'MW': 'Megafauna + Cattle Exclusion',
    'C': 'Control',
    'O': 'Cattle Only'
}

locations = {
    'OM': 'Termite Soil',
    'UT': 'Under Tree',
    'OS': 'Open Soil'
}

blocks = {
    'C': 'Central',
    'N': 'Northern',
    'S': 'Southern'
}

FLUXBOT_VOLUME = 2758   # cm^3 of fluxbot chamber
FLUXBOT_AREA = 145.5    # cm^2 area of fluxbot bottom

# Gas Law Constants:
R_m = 8.314472  # Ideal gas constant for mass [(m^3*Pa)/(K*mol)]
rho_a = 1.2041  # Density of air in kg/m^3
g_per_mol_CO2 = 44.009  # g/mol


def e_star(T, units='kPa'):
    """ Calculates saturation vapor pressure as a function of air temp.


            P = 0.61078 * exp((17.27 * T)/(T+ 273.3))

    Parameters
    ==========

            T: Air temperature, in deg-C


    Returns
    =======

            e_star: Saturation vapor pressure in kPa.

    """
    # from math import exp
    e_s = 0.61078 * np.exp((17.27 * T)/(T + 273.3))
    if units == 'kPa':
        return e_s
    elif units == 'Pa':
        return e_s*1000
    elif units == 'hPa':
        return e_s*10
    else:
        raise(ValueError, "Unknown units provided {units}".format(units=units))


def order_n_fit(x, y, n=1):
    '''Fits a polynomial fit of order n without an intercept to the data


    Parameters
    ==========

            x: np.array of time values
            y: np.array of C^prime values (CO_2(t) - CO_2(t_0))

    Returns
    =======

            p_fit: array
            - paramters for fit

            p_error: array
             standard error for the parameter estimates

            s_sq: float
                    Coefficient of determination

            f: array
                    predicted values for fit

    '''
    # Use n to determine which form of equation to fit:
    if n == 1:
        # create fitting function of form y = m*x (no b)
        def fitfunc(params, x): return params[0] * x
    elif n == 2:
        def fitfunc(params, x): return params[0] * x + params[1] * x**2
    elif n == 3:
        def fitfunc(params, x): return params[0] * \
            x + params[1] * x**2 + params[2] * x**3
    elif n > 3:
        raise(ValueError, "n must be in the range [1,3]")

    # Error function is the difference between prediction and observation:
    # create error function for least squares fit
    def errfunc(p, x, y): return fitfunc(p, x) - y

    # Initial guesses for the fitting parameters (using 0.5 for now):
    param_initial = [0.5 for i in range(n)]
    # bundle initial values in initial parameters
    init_p = np.array(tuple(param_initial))
    n_params = np.size(init_p)

    # calculate best fitting parameters (i.e. m) using the error function
    p_fit, pcov, infodict, errmsg, success = optimize.leastsq(
        errfunc, init_p.copy(), args=(x, y), full_output=1
    )  # ERROR: Unliikely, but this could be a problem? Really hope not.

    # Get fit function and set average x, y, and number of observations.
    f = fitfunc(p_fit, x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    n_obs = np.size(x)

    # sum of squared deviations from the mean for X
    ssxx = np.sum(x*x) - n_obs*x_mean*x_mean
    sse = np.sum((f - y)**2)  # residual sum of squares of errors
    # variance of y(x), predicted values of y at x versus actual values of y at x (using sse)
    syx2 = sse/(n_obs-n_params)
    sm2 = syx2 / ssxx
    ssyy = np.sum((y - y_mean)**2)

    # Determine r2 value for this fit.
    r2 = (ssyy - sse)/ssyy
    # Determine upper and lower confidence intervals
    CI_95 = 2*(math.sqrt(sm2))
    upperCI = p_fit + CI_95
    lowerCI = p_fit - CI_95

    # Determine the standard error of the estimates:
    p_error = np.absolute(pcov * syx2)**0.5

    return p_fit, p_error, ssyy, f, r2, upperCI, lowerCI


class Fluxbot():

    def __init__(self,
                 volume=FLUXBOT_VOLUME,
                 area=FLUXBOT_AREA,
                 sample_interval='60min',
                 event_type='default',
                 qaqc_max='0011',
                 smoothing_interval=20,
                 event_numbers=None,
                 filename=None,
                 output_tag='default',
                 do_avgP=False  # Uses average pressure calculated across bots.
                 ):
        """


        TODO: Add Docstrings.


                event_numbers is an optional list of events to load. This
                speeds up the loading of a fluxbot and allows for inspection of
                individual or groups of events rather than the entire dataset.

        """
        self.qaqc_max = int(qaqc_max)
        # This is our rounding interval for fluxbots.
        self.sample_interval = sample_interval
        self.event_type = event_type
        self.filename = filename
        # What interval to use when smoothing raw data for regression. Default is 5 seconds:
        self.smoothing_interval = smoothing_interval
        self._parse_filename(filename)
        # if volume == FLUXBOT_VOLUME:
        # 	print('Using default fluxbot volume {VOL}'.format(
        # 		VOL=FLUXBOT_VOLUME))
        self.volume = volume
        # if area == FLUXBOT_AREA:
        # 	print('Using default fluxbot area {AREA}'.format(
        # 		AREA=FLUXBOT_AREA))
        self.area = area
        self.df = self.read_data(filename=filename, avg_Pressure=do_avgP)
        self.events, self.bad_events = self._get_events(
            event_numbers=event_numbers)
        self.output_filename = self.filename.split('.')[0].replace(
            " ", "_") + '_output' + '_' + str(
                self.smoothing_interval) + '_' + output_tag + '.csv'

    def __repr__(self):
        return '<{0}.{1} object at {2}'.format(
            self.title, type(self).__name__, hex(id(self)))

    def _get_events(self, event_numbers=None):

        events = []
        bad_events = []
        if event_numbers:
            flux_event_numbers = event_numbers
        else:
            flux_event_numbers = self.df['Flux_Event'].unique()
        for event_number in flux_event_numbers:
            """ Here is where we create our events. They need to match
            the event type that we subclass below. Options are:

            SimpleFluxEvent (simple slope of T=0 to T=300),
            CIRASFluxEvent (quadratic),
            and LinearFluxEvent.

            The default is to use the Event class, which requires us to pass a fit_func
            to the event in order to do calculations. The default fit_func is order_2_fit

            """
            if self.event_type == 'default':
                event = Event(
                    self.df, number=event_number,
                    smoothing_interval=self.smoothing_interval
                )
            else:
                raise(ValueError, "event_type: {} is not a valid argument".format(
                    self.event_type))
            if event.is_valid(self.qaqc_max):
                events.append(event)
            else:
                bad_events.append(event)
        return events, bad_events

    def _parse_filename(self,file):

        # 1. Split the file between the directory and file name.
        [dir_name, file_name] = file.split('/')[-2:]

        # 2. The block is the first letter of the filename
        self.block = blocks[file_name[0]]

        # 3. Split the file_name to get the treatment.
        [treatment, remainder] = file_name.split('_')
        # Treament is the characters after Block
        self.treatment = treatment[1:]

        # 4. Split the remainder to get the location.
        [location, ext] = remainder.split('.')

        # 5. Replicate is the last character of the location
        self.replicate = location[-1]

        self.location = locations[location[0:-1]]

        self.title = "{block} {treatment} Plot, {location} Replicate {replicate}".format(
            block=self.block,
            treatment=self.treatment,
            location=self.location,
            replicate=self.replicate)

    def get_event(self, number=None):
        """
        Returns an event corresponding to the event number

        """
        for event in self.events:
            if event.number == number:
                return event

    def read_data(self, filename=None, avg_Pressure=False):
        """ Add docstrings later """
        df = pd.read_csv(
            filename,
            infer_datetime_format=True
        )

        # before sorting the data by timestamp, because the data were concatenated in random order, remove all the second-row
        # duplicates of column Unix.Epoch.Time. These occurred at the start of each flux event because a measurement was taken
        # at the same time right before the box closed, and right after.  Here, we're deleting the first row at
        # the start of each flux event.
        df = df.drop_duplicates(subset='Unix.Epoch.Time', keep='first')

        # wahoo! now sort.
        df = df.sort_values(by='Timestamp')
        df = df.set_index(pd.DatetimeIndex(df['Timestamp']))
        df.sort_index(inplace=True)

        # Find Events in this fluxbot's dataframe.
        # Step 1. Find the difference of ActuatorState and assign to new column
        df['StateChange'] = df.ActuatorState.diff()
        df.loc[df['StateChange'] == -0.9, 'StateChange'] = 0
        df['StateChange'] = df['StateChange'].fillna(0)
        df['Obs_Interval'] = round(
            (df['StateChange'].cumsum()/0.9)).astype('int')
        df['Flux_Event'] = 0
        df['Flux_Event'] = np.where(
            df.ActuatorState == 0, df['Obs_Interval'].values, df['Flux_Event'].values
        )

        # Smooth the data we will be using for analysis:
        # df['Raw.CO2.PPM'] = self.smooth(
        # df['Raw.CO2.PPM'],
        # smoothing_interval=self.smoothing_interval)
        if avg_Pressure is True:
            df['Pressure'] = df['avgP_perhour']

        # Unit Conversions:
        # These have been checked against the CIRAS calculations and are known to work.

        df['CO2_conc'] = self.gas_ppm_to_conc(
            df['Raw.CO2.PPM'], df['Pressure'],  df['Temp'], df['Humidity'])
        df['CO2_mass'] = self.gas_conc_to_mass(df['CO2_conc'], self.volume)

        return df

    @staticmethod
    def gas_conc_to_mass(gas_conc, volume=FLUXBOT_VOLUME):
        """
        Converts gas density in kg/m^3 to units of mass, kg.

        Parameters
        ==========

                gas_conc: panda.Series of gas concentration in units of kg/m3
                volume: float. Volume of the fluxbot chamber in cubic cm.

        Returns
        =======

                gas_mass: panda.Series of gas mass [kg]

        """
        return gas_conc * (volume * 0.000001)

    def _air_density(P, T):
        """ Calculates dry air density in kg/m^3.


        Parameters
        ==========

                P: Pressure, in hPa
                T: Temperature, in deg-C

        Returns
        =======

                rho_a, air density in kg/m^3
        """
        # Gas Law Constants:
        R_specific = 287.058  # Specific gas constant for dry air [J/kg*K]
        # Determine the density of air (assume dry air for now):
        rho_a = (P*100)/(R_specific*(T+273.15))  # Density of air in kg/m^3
        return rho_a

    def _air_mol_volume(rho_a, volume=FLUXBOT_VOLUME):
        """ Determines the moles of dry air in a given volume


        Parameters
        ==========

                rho_a: Air density, kg/m^3
                volume: Volume of chamber, cm^3


        Returns
        =======

                mol_air: Moles of air in volume.

        """
        mol_mass = 28.9628  # molar mass of dry air g/mol
        mol_density = (rho_a*1000) / mol_mass  # Molar density in mol/m^3
        # Return moles of air in volume
        return mol_density * (volume * 0.000001)

    def generate_output(self, valid_only=True):
        df_list = []
        if valid_only:
            print(
                "Only generating output for valid events (QAQC < {:04d})".format(
                    self.qaqc_max)
            )
            events = self.events
        else:
            print("Generating output for events and bad events")
            events = self.events + self.bad_events

        for event in events:
            if valid_only:
                if event.is_valid(self.qaqc_max):
                    df = event.output()
                    df['datafile'] = self.filename
                    df['chamber_volume_cm3'] = self.volume
                    df['chamber_area_cm2'] = self.area
                    df['treatment'] = self.treatment
                    df['block'] = self.block
                    df['location'] = self.location
                    df['replicate'] = self.replicate
                    df_list.append(df)
            else:
                df = event.output()
                df['datafile'] = self.filename
                df['chamber_volume_cm3'] = self.volume
                df['chamber_area_cm2'] = self.area
                df['treatment'] = self.treatment
                df['block'] = self.block
                df['location'] = self.location
                df['replicate'] = self.replicate
                df_list.append(df)
        self.output = pd.concat(df_list)

    def write(self):
        self.output.to_csv(self.output_filename, index=False)

    def _saturation_density(self, pws, pa):
        """ Returns the saturation density (kg-water/kg-air)
        given the saturation vapor pressure (pws, Pascals) and air
        pressure (pa, Pasacals)
        """
        return 0.62198 * pws / (pa - pws)

    def _humidity_ratio(self, P, T, relHum):
        """
        Determins a humidity ratio for a given pressure, temperature, and
        relative Humidity

        Parameters
        ==========

                P: panda.Series of air pressure in millibars or hPa
                T: panda.Series of air temperature in deg-C
                RH: panda.Series of relative humidity in percent [0-100]

        Returns
        =======

            x: panda.Series of humidity ration (kg-H20/kg-air)

        """
        e_s = T.apply(e_star, units='hPa')
        x_s = self._saturation_density(e_s, P)
        x = x_s * relHum / 100 # Divide by 100 to convert to fraction
        return x

    def gas_ppm_to_conc(self, gas_ppm, P, T, RelHum, g_per_mol=44.009):
        """
        Converts a gas in ppm to concentration in kg/m^3

        Parameters
        ==========

                gas_ppm: panda.Series of gas concentration in parts per million
                P: panda.Series of air pressure in millibars or hPa.
                T: panda.Series of air temperature in deg-C
                RelHum: panda.Series of relative humidity in % [0-100]
                g_per_mol: float value of the molar mass of the gas. Default is CO2 (44.009 g/mol)

        Returns
        =======
                dens: panda.Series of gas concentration in kg-CO2/m3 of air.

        """
        # Gas Law Constants:
        R_m = 8.314472  # Ideal gas constant for mass [(m^3*Pa)/(K*mol)]
        R_specific = 287.058  # Specific gas constant for dry air [J/kg*K]
        mol_mass = 28.9628  # molar mass of dry air g/mol

        # Determine the density of air, including humidity

        # Determine the density of air (assume dry air for now):
        rho_d = (P*100)/(R_specific*(T+273.15))  # Density of air in kg/m^3
        x = self._humidity_ratio(P, T, RelHum)
        rho_a = rho_d * (1 + x) / (1 + 1.609 * x)
        # Convert CO2 in ppm to CO2 in mol/m^3:
        mol_kg = 1 / (mol_mass/1000)  # Mass of air in mol/kg
        mol_m3 = mol_kg * rho_a  # mol of air per m3

        # Convert Pressure to Pa and Temperature to Kelvin.
        # P = P * 100 # hPa/millibars to Pa
        # T = T + 273.15  # Temperature in K.\

        # Use ideal gas law to determine molar density of air in mol/m^3:
        # mol_m3 = P/((R_m*T)*rho_a)

        # Convert CO2 in ppm to molar density (mol-CO2/m^3-air)
        mol_gas_m3 = mol_m3 * gas_ppm/1000000

        # Convert molar density of CO2 into CO2 density in kg/m^3:
        # CO2 concentration in kg/m^3
        gas_conc = mol_gas_m3 * (g_per_mol / 1000)

        return gas_conc  # [kg/m^3]

    @staticmethod
    def gas_ppm_to_conc_CIRAS(gas_ppm, P, T, g_per_mol=44.009):
        """
        Alternative method for converting gas in ppm to concentration in kg/m^3

        This equation is used by the CIRAS system (Eq 1., page 242-244 of CIRAS manual, v. 2018).
        It differs from gas_ppm_to_conc in that it uses the molar volume of an ideal gas at STP to
        determine molar density of gas rather than calculating the molar density of air
        directly.

        Parameters
        ==========

                gas_ppm: panda.Series of gas concentration in parts per million
                P: panda.Series of air pressure in millibars or hPa.
                T: panda.Series of air temperature in deg-C
                g_per_mol: float value of the molar mass of the gas. Default is CO2 (44.009 g/mol)

        Returns
        =======
                gas_dens: panda.Series of gas concentration in kg-gas/m3 of air


        """
        # constants:
        molarV = 22.4  # L, molar volume of an ideal gas at STP
        atm = 1013  # number of hPa per atm of pressure

        gas_conc_CIRAS = (gas_ppm) * (g_per_mol/molarV) * \
            (273.15/(T+273.15)) * (P/atm) * (1e-6)

        return gas_conc_CIRAS  # [kg/m^3]

    @staticmethod
    def plot_interval(data, title, variable='Filter.CO2.PPM'):
        """ Plots an observation interval with chosen variable

                default variable is Filter.CO2.PPM

        """
        title = title + ", {variable}".format(variable=variable)
        plt_data = [go.Scatter(x=data.index, y=data[variable])]
        layout = go.Layout(
            title=title,
            xaxis=dict(
                title='Time',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                ),
                rangeslider=dict(
                    visible=True
                )
            ),
            yaxis=dict(
                title='CO$_2$ [ppm]',
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                )
            )
        )
        fig = go.Figure(data=plt_data, layout=layout)
        plotly.offline.iplot(fig, filename=title)
        return fig


"""
Events need to calculate a flux.
We calculate a flux by fitting a curve to the data.

1. We estiamte the parameters of the curve

2. We use those parameters to determine the flux

params = Event.fit_data()

Event.calculate_flux(params)


"""


class Event():

    def __init__(self, df, number=None, smoothing_interval=20):
        """
        Must document this!

        fit_func = Fitting function to use between time and CO2 conc for flux calculation

        default fit func is the 2nd order fit.

        """
        self.raw_df = df.loc[(df.Obs_Interval == int(number))]
        self.number = number
        self.smoothing_interval = smoothing_interval
        self.smoothed_variables = ['Raw.CO2.PPM',
                                   'CO2_mass', 'CO2_conc', 'Humidity']
        self.data = Event.get_event(
            self.raw_df,
            number=self.number,
            smoothing_interval=self.smoothing_interval,
            smoothed_variables=self.smoothed_variables
        )
        self.qaqc = Event.QAQC(self.data)
        self.ambient = self.get_ambient(variable='CO2_mass')
        self.ambient_humidity = self.get_ambient(variable='Humidity')
        self.ambient_ppm = self.get_ambient(variable='Filter.CO2.PPM')
        try:
            # Raw timestamp that the event started
            self.timestamp_raw = min(self.data.index)
            # Because we usually start our measurements just before the hour (8:55-9am),
            # we also create a timestamp rounded to the nearest 30 minutes. This will help
            # with plotting our data, etc...
            self.timestamp_30_min = min(self.data.index.round(
                '30min'))  # Rounded timestamp to nearest 30 min
        except:
            self.timestamp_raw = np.nan
            self.timestamp_30_min = np.nan

        # Add the fits and parameters to this object:
        self.determine_fits()
        # Assign the change in CO2. Default is to use 2nd order slope estimate:
        self.assign_change_in_CO2()

    @staticmethod
    def smooth(variable, smoothing_interval=20):
        # smoothing mathematically using a rolling mean smoothing function
        # have tried this with windows of up to 30; we'll have to play
        # around to see what makes the most sense for our purposes.
        return variable.rolling(window=smoothing_interval, center=False).mean()

    def is_valid(self, limit):
        try:
            limit = int(limit)
        except:
            raise(ValueError, "limit: {} cannot be coerced to an integer".format(limit))
        if self.qaqc >= int(limit):
            return False
        else:
            return True

    def get_final(self, variable='Humidity', method='standard'):
        """ Get the pre-event value of a variable. If not ambient data is part of
        the data frame, return the value for the first record.

        Parameters
        ==========
                df: pandas.DataFrame of event data, including ambient if available.
                variable: string. column of event dataframe. default is CO2_mass

        Returns
        =======
                ambient in units of variable

        """
        if method == 'standard':
            try:
                final = self.data[variable][-1]
                return final
            except:
                return np.nan

    def get_ambient(self, variable='CO2_mass', method='standard'):
        """ Get the pre-event value of a variable. If not ambient data is part of
        the data frame, return the value for the first record.

        Parameters
        ==========
                df: pandas.DataFrame of event data, including ambient if available.
                variable: string. column of event dataframe. default is CO2_mass

        Returns
        =======
                ambient in units of variable

        """
        if method == 'standard':
            try:
                ambient_data = self.raw_df.loc[
                    (self.raw_df.Obs_Interval == int(self.number)) & (
                        self.raw_df.Flux_Event != int(self.number))]

                if len(ambient_data) == 0:
                    # If we only have event data, then the first observation is taken as ambient.
                    ambient = self.raw_df.iloc[0][variable]
                else:
                    ambient = ambient_data[variable][-6:-3].mean()

                return ambient
            except:
                return np.nan
        else:
            raise NotImplementedError('Method {method} is not implemented'.format(
                method=method))

    @staticmethod
    def QAQC(data, dP_max=10, dT_max=2.5, min_obs=270, dCO2_min=10, CO2_max=3000):
        """ Checks Flux Event data to see if it passes QAQC criteria


                0 0 0 0 0 0

                6 5 4 3 2 1

                1 - if dP > dP_max (default dP_max = 10 hPa)

                2 - if dT > dT_max (default dT_max = 2.5 deg-C)

                3 - if CO2 > CO2_max (default CO2_max = 3000 ppm)

                4 - if n_obs < min_obs (default min_obs = 270)

                5 - if dCO2 < dCO2_min (default dCO2_min = 10 ppm)

                6 - if sCO2 < 0 (sCO2 is difference between last CO2 ppm value and first value)

                7 - if mCO2 is True

                QAQC FLAGS:

                        1. dP (change in pressure) > 5hPa (1)
                        2. dT (change in temp) > 2 deg-C  (10)
                        3. CO2 > 3000ppm (100)
                        4. n_obs < 270 (1000)
                        5. dCO2 > 10ppm (10000)
                        6. sCO2 < 0 (100000)
                        7. mCO2 is True (1000000)

                returns qaqc, which is sum of flags.

        """
        if len(data) == 0:
            return 1000
        else:
            qaqc = 0
            dP = data['Pressure'].max()-data['Pressure'].min()
            dT = data['Temp'].max()-data['Temp'].min()
            CO2 = data['Raw.CO2.PPM'].max()
            dCO2 = data['Raw.CO2.PPM'].max() - data['Raw.CO2.PPM'].min()
            sCO2 = data['Raw.CO2.PPM'][-1] - data['Raw.CO2.PPM'][0]
            mCO2 = (
                data['Raw.CO2.PPM'][-1] > data['Raw.CO2.PPM'].mean()) & (
                data['Raw.CO2.PPM'][0] < data['Raw.CO2.PPM'].mean())
            n_obs = len(data)
            if dP > dP_max:
                qaqc = qaqc + 1
            if dT > dT_max:
                qaqc = qaqc + 10
            if CO2 > CO2_max:
                qaqc = qaqc + 100
            if n_obs < min_obs:
                qaqc = qaqc + 1000
            if dCO2 < dCO2_min:
                qaqc = qaqc + 10000
            if sCO2 < 0:
                qaqc = qaqc + 100000
            if mCO2 is False:
                qaqc = qaqc + 1000000
            return qaqc

    def determine_fits(self, variable='CO2_mass'):
        """ This function will estimate the parameters of a fit that depend on
        the order of our fitting function and save the output into the event object.

        """
        try:
            T = pd.Series(
                (self.data.index - self.data.index.min()).total_seconds())
            T.index = self.data.index
            C = self.data[variable] - self.ambient

            self.betas = []
            self.errors = []
            self.s_sqs = []
            self.f_fit = []
            self.r2_fit = []
            self.upperCIs = []
            self.lowerCIs = []
            for i in range(1, 4):
                beta_fit, error_fit, s_sq_fit, f_fit, r2_fit, upperCI, lowerCI = order_n_fit(x=T, y=C, n=i)  # NOQA
                self.betas.append(beta_fit)
                self.upperCIs.append(upperCI)
                self.lowerCIs.append(lowerCI)
                self.s_sqs.append(s_sq_fit)
                self.f_fit.append(f_fit)
                self.r2_fit.append(r2_fit)
        except:
            self.betas = np.nan
            self.upperCIs = np.nan
            self.lowerCIs = np.nan
            self.s_sqs = np.nan
            self.f_fit = np.nan
            self.r2_fit = np.nan

    def assign_change_in_CO2(self, method='max', b=None):

        try:
            # Default is to just use the beta_0 from the 2nd order polynomial.
            if method is 'max':
                # Use the maximu66m value of beta_0 from the linear and 2nd order fits:
                betas = [self.betas[0][0], self.betas[1][0]]
                beta = max(betas)  # beta in kg/s
                idx_max = max(range(len(betas)), key=betas.__getitem__)
                upperCI = self.upperCIs[idx_max][0]  # upper estimate of kg/s
                lowerCI = self.lowerCIs[idx_max][0]  # lower estimate of kg/s
            else:
                raise(ValueError, 'method must be "max"')
        except:  # If there's a problem, then just put in None
            upperCI = np.nan
            lowerCI = np.nan
            beta = np.nan

        self.beta = beta  # in kg/s
        self.upperCI = upperCI  # in kg/s (this is our upper estimate of beta)
        self.lowerCI = lowerCI  # in kg/s (this is our lower estimate of beta)

        # NOTE: We still care about the quality of this fit (i.e. R2)
        # WARNING: We assume that we only use the first parameter of our fit.
        try:
            self.duration = (self.data.index.max()
                             - self.data.index.min()).seconds
        except:
            self.duration = np.nan
        if self.beta and self.duration:
            self.change_in_kg = self.beta * self.duration  # Units of kg
        else:
            self.change_in_kg = np.nan
        if self.upperCI:
            self.change_in_kg_max = self.upperCI * \
                self.duration  # Should be in units of kg
        else:
            self.change_in_kg_max = np.nan
        if self.lowerCI:
            self.change_in_kg_min = self.lowerCI * \
                self.duration  # Should be in units of kg
        else:
            self.change_in_kg_min = np.nan

    def calculate_flux(self, area=FLUXBOT_AREA, units='kg/m2/s'):
        """

        area should be in cm2.

        Options for units are:

                'kg/m2/s'
                'g/m2/s'
                'g/cm2/s'
                'umol/m2/s'

        """

        # Convert to mass units
        try:
            change_in_g = self.change_in_kg * 1000
            change_in_g_max = self.change_in_kg_max * 1000
            change_in_g_min = self.change_in_kg_min * 1000
            # change_in_g_error = self.change_in_kg_error * 1000
            change_in_mol = change_in_g / g_per_mol_CO2
            change_in_mol_max = change_in_g_max / g_per_mol_CO2
            change_in_mol_min = change_in_g_min / g_per_mol_CO2
            # change_in_mol_error = change_in_g_error / g_per_mol_CO2
        except:
            change_in_g = np.nan
            change_in_g_max = np.nan
            change_in_g_min = np.nan
            # change_in_g_error = np.nan
            change_in_mol = np.nan
            change_in_mol_max = np.nan
            change_in_mol_min = np.nan
            # change_in_mol_error = np.nan

        # Determine areas:
        area_cm2 = area
        area_m2 = (area * 0.0001)

        try:
            if units == 'kg/m2/s':
                flux = self.change_in_kg / self.duration / area_m2  # kg / m2 /sec
                flux_max = self.change_in_kg_max / self.duration / area_m2  # kg / m2 /sec
                flux_min = self.change_in_kg_min / self.duration / area_m2  # kg / m2 /sec
                # error = self.change_in_kg_error / self.duration / area_m2
            elif units == 'g/m2/s':
                flux = change_in_g / self.duration / area_m2
                flux_max = change_in_g_max / self.duration / area_m2
                flux_min = change_in_g_min / self.duration / area_m2
                # error = change_in_g_error / self.duration / area_m2
            elif units == 'g/cm2/s':
                flux = change_in_g / self.duration / area_cm2
                flux_max = change_in_g_max / self.duration / area_cm2
                flux_min = change_in_g_min / self.duration / area_cm2
                # error = change_in_g_error / self.duration / area_cm2
            elif units == 'umol/m2/s':
                flux = (change_in_mol*1e6) / self.duration / area_m2
                flux_max = (change_in_mol_max*1e6) / self.duration / area_m2
                flux_min = (change_in_mol_min*1e6) / self.duration / area_m2
                # error = (change_in_mol_error*1e6) / self.duration / area_m2
            else:
                raise(ValueError, 'units: {} is not a valid argument'.format(units))
            return flux, flux_max, flux_min
        except:
            return np.nan, np.nan, np.nan

    def output(self):
        try:
            _1st_order_beta_0 = self.betas[0][0]
        except:
            _1st_order_beta_0 = np.nan
        try:
            _1st_order_beta_0_max = self.upperCIs[0][0]
        except:
            _1st_order_beta_0_max = np.nan
        try:
            _1st_order_beta_0_min = self.lowerCIs[0][0]
        except:
            _1st_order_beta_0_min = np.nan
        try:
            _1st_order_r_sq = self.r2_fit[0]
        except:
            _1st_order_r_sq = np.nan
        try:
            _2nd_order_beta_0 = self.betas[1][0]
        except:
            _2nd_order_beta_0 = np.nan
        try:
            _2nd_order_beta_0_max = self.upperCIs[1][0]
        except:
            _2nd_order_beta_0_max = np.nan
        try:
            _2nd_order_beta_0_min = self.lowerCIs[1][0]
        except:
            _2nd_order_beta_0_min = np.nan
        try:
            _2nd_order_beta_1 = self.betas[1][1]
        except:
            _2nd_order_beta_1 = np.nan
        try:
            _2nd_order_beta_1_max = self.upperCIs[1][1]
        except:
            _2nd_order_beta_1_max = np.nan
        try:
            _2nd_order_beta_1_min = self.lowerCIs[1][1]
        except:
            _2nd_order_beta_1_min = np.nan
        try:
            _2nd_order_r_sq = self.r2_fit[1]
        except:
            _2nd_order_r_sq = np.nan
        event_data = {}
        flux, flux_max, flux_min = self.calculate_flux(units='umol/m2/s')
        if type(self.timestamp_raw) is pd.Timestamp:
            event_data['timestamp'] = self.timestamp_raw.isoformat()
        else:
            event_data['timestamp'] = pd.Timestamp.now()
        if type(self.timestamp_30_min) is pd.Timestamp:
            event_data['year'] = self.timestamp_30_min.year
            event_data['month'] = self.timestamp_30_min.month
            event_data['day'] = self.timestamp_30_min.day
            event_data['hour'] = self.timestamp_30_min.hour
        else:
            event_data['year'] = '{:5.1f}'.format(np.nan)
            event_data['month'] = '{:5.1f}'.format(np.nan)
            event_data['day'] = '{:5.1f}'.format(np.nan)
            event_data['hour'] = '{:5.1f}'.format(np.nan)
        event_data['event_number'] = self.number
        event_data['avg_temp_degC'] = '{:5.1f}'.format(self.data.Temp.mean())
        event_data['avg_pressure_hPa'] = '{:5.1f}'.format(
            self.data.Pressure.mean())
        event_data['avg_rel_humidity'] = '{:5.1f}'.format(
            self.data.Humidity.mean())
        event_data['ambient_humidity'] = '{:5.1f}'.format(
            self.ambient_humidity)
        event_data['final_humidity'] = '{:5.1f}'.format(
            self.get_final('Humidity'))
        event_data['ambient_CO2_kg'] = '{:5.2E}'.format(self.ambient)
        event_data['ambient_CO2_ppm'] = '{:5.0f}'.format(self.ambient_ppm)
        event_data['beta'] = '{:5.2E}'.format(self.beta)
        event_data['duration_sec'] = self.duration
        event_data['change_in_CO2_kg'] = '{:5.2E}'.format(self.change_in_kg)
        event_data['1st_order_beta_0'] = '{:5.2E}'.format(_1st_order_beta_0)
        event_data['1st_order_beta_0_max'] = '{:5.2E}'.format(
            _1st_order_beta_0_max)
        event_data['1st_order_beta_0_min'] = '{:5.2E}'.format(
            _1st_order_beta_0_min)
        event_data['1st_order_r_sq'] = '{:5.3f}'.format(_1st_order_r_sq)
        event_data['2nd_order_beta_0'] = '{:5.3E}'.format(_2nd_order_beta_0)
        event_data['2nd_order_beta_0_max'] = '{:5.2E}'.format(
            _2nd_order_beta_0_max)
        event_data['2nd_order_beta_0_min'] = '{:5.2E}'.format(
            _2nd_order_beta_0_min)
        event_data['2nd_order_beta_1'] = '{:5.2E}'.format(_2nd_order_beta_1)
        event_data['2nd_order_beta_1_max'] = '{:5.2E}'.format(
            _2nd_order_beta_1_max)
        event_data['2nd_order_beta_1_min'] = '{:5.2E}'.format(
            _2nd_order_beta_1_min)
        event_data['2nd_order_r_sq'] = '{:5.3f}'.format(_2nd_order_r_sq)
        event_data['flux_umol_m2_sec'] = '{:5.3f}'.format(flux)
        # event_data['flux_error_umol_m2_sec'] = '{:5.3f}'.format(error)
        event_data['flux_max_umol_m2_sec'] = '{:5.3f}'.format(flux_max)
        event_data['flux_min_umol_m2_sec'] = '{:5.3f}'.format(flux_min)
        event_data['qaqc_flags'] = '{:04d}'.format(self.qaqc)
        return pd.DataFrame(event_data, index=[0])

    def write(self):
        """ Special function to write out only the data for this event.

        We write out the processed, smoothed data to ensure consistency.

        """
        return None

    def plot_event(self, variable='Filter.CO2.PPM'):
        """ Plots a specific event number

        """
        title = 'Raw Data for Event #{num}'.format(num=self.number)
        Fluxbot.plot_interval(self.data, title, variable=variable)

    @ staticmethod
    def get_event(df, number=None, smoothing_interval=None, smoothed_variables=None):
        """ Returns a new dataframe containing a specific event
                from a fluxbot dataset.

                number = The event number from the fluxbot dataset
                cutoff = The number of observations (seconds) to remove from the start
                                 of the dataframe. Used to handle memory/initial condition
                                 effects in the CO2 data.

        """
        # Step 1. Truncate dataframe to only include Event data (get rid of ambient data)
        df = df.loc[(df.Flux_Event == int(number))]

        # Step 2. If necessary, smooth event data.
        if smoothing_interval:
            if smoothed_variables:
                for variable in smoothed_variables:
                    df[variable] = Event.smooth(
                        df[variable], smoothing_interval=smoothing_interval)
            df = df[smoothing_interval:]
        return df
        return df
        return df
        return df
        return df
