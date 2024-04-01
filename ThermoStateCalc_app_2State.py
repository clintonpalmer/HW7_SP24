#
#   Description:
#               A graphical user interface for calculating the thermodynamic properties of steam at two different states
#                   by specifying any two properties (e.g., p&T, s&h, p&v, etc) at each state.
#
#
#   Assumptions:
#               You've read the Assumptions, Instructions and Important Info below.
#               You will enter reasonable values into the calculator, otherwise you will likely get an error.
#               A negative result in State Change implies a decrease in the value from State 1 to State 2.
#
#   Instructions:
#                   1. Press start button.
#                   2. User interface will display.
#                   3. Select properties for State 1.
#                   4. Enter --specific-- values for State 1.
#                   5. press "Calculate State 1" button  (See "Important Info" below)
#                       -results for State 1 appear
#                   6. repeat for State 2
#                       -results for State 2 and State Change appear. (Note:  the results for State Change update
#                                                                            only when "Calculate State 2 / State Change" button
#                                                                            is pressed)
#
#
#   Important Info:
#                   Updating the "State Change" data is tied to the "Calculate State 2 / State Change" button.  If
#                    you open the GUI and press "Calculate State 2" button first, the GUI will close.  When you
#                    first open the GUI, you must perform the State 1 calculations before calculating State 2.
#                    After pressing the "Calculate State 1" button for the first time, the order does not matter.
#
# Chatgpt assisted with the developement of this code.


#region imports
import sys
from ThermoStateCalc_2State import Ui__frm_StateCalculator
from pyXSteam.XSteam import XSteam
from PyQt5.QtWidgets import QWidget, QApplication
from UnitConversion import UC
from scipy.optimize import fsolve
#endregion

#region class definitions
class main_window(QWidget,Ui__frm_StateCalculator):
    """
    This class represents the main window of the state calculator application.
    It inherits from QWidget and Ui__frm_StateCalculator to incorporate both
    the QWidget functionality and the UI design created with Qt Designer.

    Attributes:
        steamTable (XSteam): Instance of the XSteam class for steam property calculations.
        currentUnits (str): String representing the current unit system ('SI' or 'English').
    """
    def __init__(self):
        """
        Initializes the main window.

        This function sets up the UI, connects signals and slots, initializes the steam table,
        sets the current units to 'SI', applies the units setup, and shows the window.
        """
        super().__init__()
        self.setupUi(self)
        self.SetupSlotsAndSignals()
        self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.currentUnits='SI'
        self.setUnits()
        self.show()

    def SetupSlotsAndSignals(self):
        """
        This function sets up connections between specific user interface elements and corresponding
        functions to handle user interactions. These connections ensure that certain actions are
        performed when the user interacts with specific parts of the interface.

        Signals Connected:
        - Clicked signals of radio buttons _rdo_English and _rdo_SI are connected to the setUnits method.
        - CurrentIndexChanged signals of combo boxes _cmb_Property1 and _cmb_Property2 are connected to the setUnits method.
        - Clicked signal of push button _pb_Calculate is connected to the calculateProperties method.
        - CurrentIndexChanged signals of combo boxes _cmb_Property1_2 and _cmb_Property2_2 are connected to the setUnits method.
        - Clicked signal of push button _pb_Calculate_2 is connected to the calculatestate2 method.
        """
        self._rdo_English.clicked.connect(self.setUnits)
        self._rdo_SI.clicked.connect(self.setUnits)
        self._cmb_Property1.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate.clicked.connect(self.calculateProperties)
        self._cmb_Property1_2.currentIndexChanged.connect(self.setUnits)
        self._cmb_Property2_2.currentIndexChanged.connect(self.setUnits)
        self._pb_Calculate_2.clicked.connect(self.calculatestate2)
        pass

    def setUnits(self):
        """
        This sets the units for the selected specified properties.
        Units for the thermodynamic properties are set upon pushing calculate button.
        :return:
        """
        #set the units system based on selected radio button
        #also, determine if a units change is required
        SI=self._rdo_SI.isChecked()
        newUnits='SI' if SI else 'EN'
        UnitChange = self.currentUnits != newUnits  # compare new units to current units
        self.currentUnits = newUnits

        if SI:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_MKS)
            self.l_Units = "m"
            self.p_Units = "bar"
            self.t_Units = "C"
            self.m_Units = "kg"
            self.time_Units = "s"
            self.energy_Units = "W"
            self.u_Units = "kJ/kg"
            self.h_Units = "kJ/kg"
            self.s_Units = "kJ/kg*C"
            self.v_Units = "m^3/kg"
        else:
            self.steamTable=XSteam(XSteam.UNIT_SYSTEM_FLS)
            self.l_Units = "ft"
            self.p_Units = "psi"
            self.t_Units = "F"
            self.m_Units = "lb"
            self.time_Units = "s"
            self.energy_Units = "btu"
            self.u_Units = "btu/lb"
            self.h_Units = "btu/lb"
            self.s_Units = "btu/lb*F"
            self.v_Units = "ft^3/lb"

        #read selected Specified Properties from combo boxes
        SpecifiedProperty1 = self._cmb_Property1.currentText()
        SpecifiedProperty2 = self._cmb_Property2.currentText()
        SpecifiedProperty1_2 = self._cmb_Property1_2.currentText()
        SpecifiedProperty2_2 = self._cmb_Property2_2.currentText()
        #read numerical values for selected properties
        SP=[float(self._le_Property1.text()), float(self._le_Property2.text()),
            float(self._le_Property1_2.text()), float(self._le_Property2_2.text())]

        #set units labels and convert values if needed
        if 'Pressure' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.p_Units)
            if UnitChange:  # note that I only should convert if needed.  Not if I double click on SI or English
                SP[0]=SP[0]*UC.psi_to_bar if SI else SP[0]*UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.t_Units)
            if UnitChange:
                SP[0] = UC.F_to_C(SP[0]) if SI else UC.C_to_F(SP[0])
        elif 'Energy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.u_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.btuperlb_to_kJperkg if SI else SP[0]*UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.h_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlb_to_kJperkg if SI else SP[0] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.s_Units)
            if UnitChange:
                SP[0] = SP[0] * UC.btuperlbF_to_kJperkgC if SI else SP[0] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText(self.v_Units)
            if UnitChange:
                SP[0]=SP[0]*UC.ft3perlb_to_m3perkg if SI else SP[0]*UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty1:
            self._lbl_Property1_Units.setText("")


        if 'Pressure' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.p_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.psi_to_bar if SI else SP[1] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.t_Units)
            if UnitChange:
                SP[1] = UC.F_to_C(SP[1]) if SI else UC.C_to_F(SP[1])
        elif 'Energy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.u_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.h_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlb_to_kJperkg if SI else SP[1] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.s_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.btuperlbF_to_kJperkgC if SI else SP[1] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText(self.v_Units)
            if UnitChange:
                SP[1] = SP[1] * UC.ft3perlb_to_m3perkg if SI else SP[1] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty2:
            self._lbl_Property2_Units.setText("")

        if 'Pressure' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.p_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.psi_to_bar if SI else SP[2] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.t_Units)
            if UnitChange:
                SP[2] = UC.F_to_C(SP[2]) if SI else UC.C_to_F(SP[2])
        elif 'Energy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.u_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlb_to_kJperkg if SI else SP[2] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.h_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlb_to_kJperkg if SI else SP[2] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.s_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.btuperlbF_to_kJperkgC if SI else SP[2] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText(self.v_Units)
            if UnitChange:
                SP[2] = SP[2] * UC.ft3perlb_to_m3perkg if SI else SP[2] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty1_2:
            self._lbl_Property1_2_Units.setText("")

        if 'Pressure' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.p_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.psi_to_bar if SI else SP[3] * UC.bar_to_psi
        elif 'Temperature' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.t_Units)
            if UnitChange:
                SP[3] = UC.F_to_C(SP[3]) if SI else UC.C_to_F(SP[3])
        elif 'Energy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.u_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlb_to_kJperkg if SI else SP[3] * UC.kJperkg_to_btuperlb
        elif 'Enthalpy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.h_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlb_to_kJperkg if SI else SP[3] * UC.kJperkg_to_btuperlb
        elif 'Entropy' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.s_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.btuperlbF_to_kJperkgC if SI else SP[3] * UC.kJperkgC_to_btuperlbF
        elif 'Volume' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText(self.v_Units)
            if UnitChange:
                SP[3] = SP[3] * UC.ft3perlb_to_m3perkg if SI else SP[3] * UC.m3perkg_to_ft3perlb
        elif 'Quality' in SpecifiedProperty2_2:
            self._lbl_Property2_2_Units.setText("")

        self._le_Property1.setText("{:0.3f}".format(SP[0]))
        self._le_Property2.setText("{:0.3f}".format(SP[1]))
        self._le_Property1_2.setText("{:0.3f}".format(SP[2]))
        self._le_Property2_2.setText("{:0.3f}".format(SP[3]))
    def clamp(self, x, low, high):
        """
        This clamps a float x between a high and low limit inclusive
        :param x:
        :param low:
        :param high:
        :return:
        """
        if x<low:
            return low
        if x>high:
            return high
        return x

    def between(self, x, low, high):
        """
        Tells if x is between low and high inclusive
        :param x:
        :param low:
        :param high:
        :return:
        """
        if x>=low and x<=high:
            return True
        return False

    def getSatProps_p(self, p):
        """
        Given a pressure, calculate the saturated properties for that isobar
        :param p:
        :return:
        """
        self.tSat=self.steamTable.tsat_p(p)
        self.pSat=p
        self.vf=self.steamTable.vL_p(p)
        self.vg=self.steamTable.vV_p(p)
        self.hf=self.steamTable.hL_p(p)
        self.hg=self.steamTable.hV_p(p)
        self.uf=self.steamTable.uL_p(p)
        self.ug=self.steamTable.uV_p(p)
        self.sf=self.steamTable.sL_p(p)
        self.sg=self.steamTable.sV_p(p)
        self.vgf=self.vg-self.vf
        self.hgf=self.hg-self.hf
        self.sgf=self.sg-self.sf
        self.ugf=self.ug-self.uf

    def getSatProps_t(self, t):
        """
        Given a temperature, calculate the saturation pressure and then
        calculate all other saturated properties
        :param t:
        :return:
        """
        self.tSat=t
        self.pSat=self.steamTable.psat_t(t)
        self.getSatProps_p(self.pSat)

    def makeLabel_1Phase(self):
        """
        For State 1:  Given temperature (T) and pressure (P), this method calculates and updates various other properties
        including region, pressure, temperature, saturation temperature, internal energy, enthalpy, entropy,
        specific volume, and quality. It then sets the label to display these properties.
        The values of each property are stored in self.property_values# and are used in the difference calculation
        Returns:
        None
        """
        self._lbl_State.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        self.u = self.steamTable.u_pt(self.p, self.t)
        self.h = self.steamTable.h_pt(self.p, self.t)
        self.s = self.steamTable.s_pt(self.p, self.t)
        self.v = self.steamTable.v_pt(self.p, self.t)
        self.x = 1.0 if self.t > self.steamTable.tsat_p(self.p) else 0.0
        stProps += "\nInternal Energy = {:0.3f} ({:})".format(self.u, self.u_Units)
        stProps += "\nEnthalpy = {:0.3f} ({:})".format(self.h, self.h_Units)
        stProps += "\nEntropy = {:0.3f} ({:})".format(self.s, self.s_Units)
        stProps += "\nSpecific Volume = {:0.3f} ({:})".format(self.v, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps
        # Create a list to store the property values
        self.property_values1 = [self.p, self.t, self.tSat, self.u, self.h, self.s, self.v, self.x]
        # Store the new values
        #print(stProps)

    def makeLabel_2Phase(self):
        """
        For State 1:  Given pressure (P) and quality (x), this method calculates and updates various other properties
        including region, pressure, temperature, saturation temperature, internal energy, enthalpy, entropy,
        specific volume, and quality. It then sets the label to display these properties.
        The values of each property are stored in self.property_values# and are used in the difference calculation
        Returns:
        None
        """
        self._lbl_State.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        stProps += "\nInternal Energy = u={:0.3f} ({:})".format(self.uf + self.x * self.ugf, self.u_Units)
        stProps += "\nEnthalpy = h={:0.3f} ({:})".format(self.hf + self.x * self.hgf, self.h_Units)
        stProps += "\nEntropy = s={:0.3f} ({:})".format(self.sf + self.x * self.sgf, self.s_Units)
        stProps += "\nSpecific Volume = v={:0.5f} ({:})".format(self.vf + self.x * self.vgf, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.property_values1 = [self.p, self.t, self.tSat, self.u, self.h, self.s, self.v, self.x]
        self.stProps=stProps
       # print(stProps)


    def makeLabel_1Phase2(self):
        """
        For State 2:  Given temperature (T) and pressure (P), this method calculates and updates various other properties
        including region, pressure, temperature, saturation temperature, internal energy, enthalpy, entropy,
        specific volume, and quality. It then sets the label to display these properties.
        The values of each property are stored in self.property_values# and are used in the difference calculation
        Returns:
        None
        """
        self._lbl_State_2.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        self.u = self.steamTable.u_pt(self.p, self.t)
        self.h = self.steamTable.h_pt(self.p, self.t)
        self.s = self.steamTable.s_pt(self.p, self.t)
        self.v = self.steamTable.v_pt(self.p, self.t)
        self.x = 1.0 if self.t > self.steamTable.tsat_p(self.p) else 0.0
        stProps += "\nInternal Energy = {:0.3f} ({:})".format(self.u, self.u_Units)
        stProps += "\nEnthalpy = {:0.3f} ({:})".format(self.h, self.h_Units)
        stProps += "\nEntropy = {:0.3f} ({:})".format(self.s, self.s_Units)
        stProps += "\nSpecific Volume = {:0.3f} ({:})".format(self.v, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps=stProps
        self.property_values2 = [self.p, self.t, self.tSat, self.u, self.h, self.s, self.v, self.x]
        self.calculate_difference()
        #print(self.property_values1)

        #print(stProps)


    def makeLabel_2Phase2(self):
        """
        For State 2:  Given pressure (P) and quality (x), this method calculates and updates various other properties
        including region, pressure, temperature, saturation temperature, internal energy, enthalpy, entropy,
        specific volume, and quality. It then sets the label to display these properties.
        The values of each property are stored in self.property_values# and are used in the difference calculation
        Returns:
        None
        """
        self._lbl_State_2.setText("Region = {:}".format(self.region))
        stProps = "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
        stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat,
                                                                               self.t_Units)
        stProps += "\nInternal Energy = u={:0.3f} ({:})".format(self.uf + self.x * self.ugf, self.u_Units)
        stProps += "\nEnthalpy = h={:0.3f} ({:})".format(self.hf + self.x * self.hgf, self.h_Units)
        stProps += "\nEntropy = s={:0.3f} ({:})".format(self.sf + self.x * self.sgf, self.s_Units)
        stProps += "\nSpecific Volume = v={:0.5f} ({:})".format(self.vf + self.x * self.vgf, self.v_Units)
        stProps += "\nQuality = {:0.3f}".format(self.x)
        self.stProps = stProps
        self.property_values2 = [self.p, self.t, self.tSat, self.u, self.h, self.s, self.v, self.x]
        self.calculate_difference()
       # print(stProps)

    def calculate_difference(self):
        """
        Calculate the difference between the "State 1" and "State 2" properties.
        :return: The difference for each property
        """
        difference = "\nPressure = {:0.3f} ({:})".format(self.property_values2[0] - self.property_values1[0],
                                                                    self.p_Units)
        difference += "\nTemperature = {:0.3f} ({:})".format(
            self.property_values2[1] - self.property_values1[1], self.t_Units)
        difference += "\nInternal Energy = {:0.3f} ({:})".format(
            self.property_values2[3] - self.property_values1[3], self.u_Units)
        difference += "\nEnthalpy = {:0.3f} ({:})".format(
            self.property_values2[4] - self.property_values1[4], self.h_Units)
        difference += "\nEntropy = {:0.3f} ({:})".format(self.property_values2[5] - self.property_values1[5],
                                                                    self.s_Units)
        difference += "\nSpecific Volume = {:0.3f} ({:})".format(
            self.property_values2[6] - self.property_values1[6], self.v_Units)
        difference += "\nQuality = {:0.3f}".format(self.property_values2[7] - self.property_values1[7])

        self._lbl_StateChangeProperties.setText(difference)
    def calculateProperties(self):
        """
        For State 1:  Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above
        :return: nothing
        """
        # Step 1: read which properties are being specified from the combo boxes
        SP=[self._cmb_Property1.currentText()[-2:-1], self._cmb_Property2.currentText()[-2:-1], self._cmb_Property1_2.currentText()[-2:-1], self._cmb_Property2_2.currentText()[-2:-1]]
        if SP[0]==SP[1]:
            self._lbl_Warning.setText("Warning:  You cannot specify the same property twice.")
        else:
            self._lbl_Warning.setText("")
        SP[0]=SP[0].lower()
        SP[1]=SP[1].lower()

        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if SP[0]=='p' or SP[1]=='p':
            oFlipped= SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 1:  pt or tp
            if SP1=='t':
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2
                self.t= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.tSat=round(self.tSat,3)  # I will compare at 3 three decimal places
                # compare T to TSat
                if self.t<self.tSat or self.t>self.tSat:
                    self.region = "sub-cooled liquid" if self.t<self.tSat else "super-heated vapor"
                    self.makeLabel_1Phase()
                else:  #this is ambiguous since at saturated temperature
                    self.region = "two-phase"
                    self.stProps ="Region = {:}".format(self.region)
                    self.stProps += "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
                    self.stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units, self.tSat, self.t_Units)
                    self.stProps += "\nInternal Energy = uf={:0.3f}, ug={:0.3f} ({:})".format(self.uf, self.ug, self.u_Units)
                    self.stProps += "\nEnthalpy = hf={:0.3f}, hg={:0.3f} ({:})".format(self.hf,self.hg, self.h_Units)
                    self.stProps += "\nEntropy = sf={:0.3f}, sg={:0.3f} ({:})".format(self.sf, self.sg, self.s_Units)
                    self.stProps += "\nSpecific Volume = vf={:0.3f}, vg={:0.3f} ({:})".format(self.vf, self.vg, self.v_Units)
                    self.stProps += "\nQuality = unknown"
        #case 2: pv or vp
            elif SP1=='v':
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2
                self.v= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find T
                    dt=1.0 if self.v>self.vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 3 pu or up
            elif SP1=='u':
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2;
                self.u= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u<self.uf or self.u>self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find T
                    dt=1.0 if self.u>self.ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 4 ph or hp
            elif SP1 == 'h':
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2;
                self.h= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.t=self.steamTable.t_ph(self.p, self.h)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.t=self.tSat
                    self.makeLabel_2Phase()
        #case 5 ps or sp
            elif SP1 == 's':
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2;
                self.s= f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.t=self.steamTable.t_ps(self.p, self.s)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
        #case 6 px or xp
            elif SP1 == 'x':
                self.region="two-phase"
                f1=float(self._le_Property1.text())
                f2=float(self._le_Property2.text())
                self.p= f1 if not oFlipped else f2
                self.x= f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.t=self.tSat
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0]=='t' or SP[1]=='t':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
        #case 7:  tv or vt
            if SP1=='v':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.t = f1 if not oFlipped else f2
                self.v = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                self.vf=round(self.vf,5)
                self.vg=round(self.vg,3)
                # compare v to vf and vg
                if self.v<self.vf or self.v>self.vg:
                    self.region = "sub-cooled liquid" if self.v<self.vf else "super-heated vapor"
                    #since I can't find properties using v, I will use fsolve to find P
                    dp=-0.1 if self.v>self.vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.v-self.vf)/(self.vgf)
                    self.p=self.pSat
                    self.makeLabel_2Phase()
        #case 8:  tu or ut
            elif SP1=='u':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.t = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u<self.uf else "super-heated vapor"
                    #since I can't find properties using u, I will use fsolve to find P
                    dp=0.1 if self.u>self.ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(self.t, P)
                    self.p = fsolve(fn8, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  #two-phase
                    self.region = "two-phase"
                    #first calculate quality
                    self.x=(self.u-self.uf)/(self.ugf)
                    self.p=self.pSat
                    self.makeLabel_2Phase()
        #case 9:  th or ht
            elif SP1 == 'h':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.t = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.p=self.steamTable.p_th(self.t, self.h)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.makeLabel_2Phase()
        #case 10:  ts or st
            elif SP1 == 's':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.t = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.p=self.steamTable.p_ts(self.t, self.s)
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p=self.pSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
        #case 11:  tx or xt
            elif SP1 == 'x':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.t = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region="two-phase"
                self.getSatProps_t(self.t)
                self.p=self.pSat
                self.x = float(self._le_Property2.text())
                self.x=self.clamp(self.x,0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 12:  vh or hv
            if SP1 == 'h':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.v = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                def fn12(P):
                    """
                    Calculates the difference between specific volume based on enthalpy and pressure.

                    This function calculates the difference between specific volume based on enthalpy and pressure.
                    It determines if the state is single-phase or two-phase based on whether the enthalpy falls
                    within the saturated liquid and vapor enthalpy range. If it's two-phase, it calculates the
                    quality (x) and returns the difference between the specific volume and the saturated mixture
                    specific volume. If it's single-phase, it calculates the specific volume based on the enthalpy
                    and pressure and returns the difference between the calculated specific volume and the one
                    obtained from the steam table.

                    Parameters:
                        P (float): Pressure value.

                    Returns:
                        float: Difference between specific volume.
                    """
                    # could be single phase or two-phase, but both v&h have to match at same x
                    self.getSatProps_p(P)
                    if self.between(self.h,self.hf, self.hg):
                        self.x=(self.h-self.hf)/self.hgf
                        return self.v-(self.vf+self.x*self.vgf)
                    # could be single phase
                    return self.v-self.steamTable.v_ph(P,self.h)
                self.p=fsolve(fn12,[1.0])[0]
                self.t=self.steamTable.t_ph(self.p, self.h)
                self.getSatProps_p(self.p)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.makeLabel_2Phase()
            # case 13:  vu or uv
            elif SP1 == 'u':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.v = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn13(PT):
                    """
                    Calculates differences in specific volume and internal energy based on pressure and temperature.

                    This function calculates the differences in specific volume and internal energy based on the given
                    pressure and temperature. It determines if the state is two-phase based on whether the internal
                    energy falls within the saturated liquid and vapor internal energy range. If it's two-phase, it
                    calculates the saturation temperature, quality (x), and returns the differences between specific
                    volume and internal energy from the saturated mixture and the actual values. If it's single-phase,
                    it calculates specific volume and internal energy based on the pressure and temperature and returns
                    the differences between the calculated values and the ones obtained from the steam table.

                    Parameters:
                        PT (list): List containing pressure and temperature values.

                    Returns:
                        list: Differences in specific volume and internal energy.
                    """
                    self.getSatProps_p(PT[0])
                    if self.between(self.u,self.uf, self.ug):
                        self.t=self.tSat
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.v-(self.vf+self.x*self.vgf),0]
                    return [self.v-self.steamTable.v_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn13, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase()
            # case 14:  vs os sv
            elif SP1 == 's':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.v = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn14(PT):
                    """
                        Calculates differences in specific volume and entropy based on pressure and temperature.

                        This function calculates the differences in specific volume and entropy based on the given
                        pressure and temperature. It determines if the state is two-phase based on whether the entropy
                        falls within the saturated liquid and vapor entropy range. If it's two-phase, it calculates
                        the quality (x) and returns the differences between specific volume and entropy from the
                        saturated mixture and the actual values. If it's single-phase, it calculates specific volume
                        and entropy based on the pressure and temperature and returns the differences between the
                        calculated values and the ones obtained from the steam table.

                        Parameters:
                            PT (list): List containing pressure and temperature values.

                        Returns:
                            list: Differences in specific volume and entropy.
                        """
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.v-self.vf-self.x*self.vgf, 0.0]
                    return [self.v - self.steamTable.v_pt(PT[0], PT[1]),
                                  self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn14, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 15:  vx or xv
            elif SP1 == 'x':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn15(p):
                    """
                    Calculates the difference in specific volume based on pressure.

                    This function calculates the difference in specific volume based on the given pressure.
                    It determines if the state is two-phase based on the quality (x) and returns the difference
                    between the specific volume and the saturated mixture specific volume. If it's single-phase,
                    it calculates specific volume based on the pressure and returns the difference between the
                    calculated value and the one obtained from the steam table.

                    Parameters:
                        p (float): Pressure value.

                    Returns:
                        float: Difference in specific volume.
                    """
                    self.getSatProps_p(p)
                    return self.v -(self.vf+self.x*self.vgf)
                props=fsolve(fn15, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 16:  hu or uh
            if SP1 == 'u':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.h = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                # use fsolve to fing P&T at this v & u
                def fn16(PT):
                    """
                    Calculates differences in enthalpy and internal energy based on pressure and temperature.

                    This function calculates the differences in enthalpy and internal energy based on the given
                    pressure and temperature. It determines if the state is two-phase based on whether the internal
                    energy falls within the saturated liquid and vapor internal energy range. If it's two-phase, it
                    calculates the quality (x) and returns the differences between enthalpy and internal energy from
                    the saturated mixture and the actual values. If it's single-phase, it calculates enthalpy and
                    internal energy based on the pressure and temperature and returns the differences between the
                    calculated values and the ones obtained from the steam table.

                    Parameters:
                        PT (list): List containing pressure and temperature values.

                    Returns:
                        list: Differences in enthalpy and internal energy.
                    """
                    self.getSatProps_p(PT[0])
                    if self.between(self.u, self.uf, self.ug):
                        self.x=(self.u-self.uf)/self.ugf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.u-self.steamTable.u_pt(PT[0],PT[1])]
                props=fsolve(fn16, [1,100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t=self.tSat
                    self.makeLabel_2Phase()
            # case 17:  hs or sh
            elif SP1 == 's':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.h = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn17(PT):
                    """
                    Calculates differences in enthalpy and entropy based on pressure and temperature.

                    This function calculates the differences in enthalpy and entropy based on the given
                    pressure and temperature. It determines if the state is two-phase based on whether the entropy
                    falls within the saturated liquid and vapor entropy range. If it's two-phase, it calculates
                    the quality (x) and returns the differences between enthalpy and entropy from the saturated
                    mixture and the actual values. If it's single-phase, it calculates enthalpy and entropy
                    based on the pressure and temperature and returns the differences between the calculated
                    values and the ones obtained from the steam table.

                    Parameters:
                        PT (list): List containing pressure and temperature values.

                    Returns:
                        list: Differences in enthalpy and entropy.
                    """
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x=(self.s-self.sf)/self.sgf
                        return [self.h-self.hf-self.x*self.hgf, 0.0]
                    return [self.h-self.steamTable.h_pt(PT[0],PT[1]),self.s-self.steamTable.s_pt(PT[0],PT[1])]
                props = fsolve(fn17, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 18:  hx or xh
            elif SP1 == 'x':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn18(p):
                    """
                    Calculates the difference in enthalpy based on pressure.

                    This function calculates the difference in enthalpy based on the given pressure.
                    It determines if the state is two-phase based on the quality (x) and returns the difference
                    between enthalpy and the saturated mixture enthalpy. If it's single-phase, it calculates
                    enthalpy based on the pressure and returns the difference between the calculated value
                    and the one obtained from the steam table.

                    Parameters:
                        p (float): Pressure value.

                    Returns:
                        float: Difference in enthalpy.
                    """
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn18, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 19:  us or su
            if SP1 == 's':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.u = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                def fn19(PT):
                    """
                    Calculates differences in internal energy and entropy based on pressure and temperature.

                    This function calculates the differences in internal energy and entropy based on the given
                    pressure and temperature. It determines if the state is two-phase based on whether the entropy
                    falls within the saturated liquid and vapor entropy range. If it's two-phase, it calculates
                    the quality (x) and returns the differences between internal energy and entropy from the
                    saturated mixture and the actual values. If it's single-phase, it calculates internal energy
                    and entropy based on the pressure and temperature and returns the differences between the
                    calculated values and the ones obtained from the steam table.

                    Parameters:
                        PT (list): List containing pressure and temperature values.

                    Returns:
                        list: Differences in internal energy and entropy.
                    """
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.u - self.uf - self.x * self.ugf, 0.0]
                    return [self.u - self.steamTable.u_pt(PT[0], PT[1]),
                            self.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn19, [1, 100])
                self.p=props[0]
                self.t=props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t=self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase()
            # case 20:  ux or xu
            elif SP1 == 'x':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.u = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn20(p):
                    """
                    Calculates the difference in enthalpy based on pressure.

                    This function calculates the difference in enthalpy based on the given pressure.
                    It determines if the state is two-phase based on the quality (x) and returns the difference
                    between enthalpy and the saturated mixture enthalpy. If it's single-phase, it calculates
                    enthalpy based on the pressure and returns the difference between the calculated value
                    and the one obtained from the steam table.

                    Parameters:
                        p (float): Pressure value.

                    Returns:
                        float: Difference in enthalpy.
                    """
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn20, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2.text()),0.0,1.0)
                self.makeLabel_2Phase()
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 21:  sx or xs
            if SP1 == 'x':
                f1 = float(self._le_Property1.text())
                f2 = float(self._le_Property2.text())
                self.s = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                def fn21(p):
                    """
                    Calculates the difference in enthalpy based on pressure.

                    This function calculates the difference in enthalpy based on the given pressure.
                    It determines if the state is two-phase based on the quality (x) and returns the difference
                    between enthalpy and the saturated mixture enthalpy. If it's single-phase, it calculates
                    enthalpy based on the pressure and returns the difference between the calculated value
                    and the one obtained from the steam table.

                    Parameters:
                        p (float): Pressure value.

                    Returns:
                        float: Difference in enthalpy.
                    """
                    self.getSatProps_p(p)
                    return self.h -(self.hf+self.x*self.hgf)
                props=fsolve(fn21, [1])
                self.p=props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2.text()),0.0,1.0)
                self.makeLabel_2Phase()
        self._lbl_StateProperties.setText(self.stProps)


    def calculatestate2(self):
        """
        For State 2:  Calculates the thermodynamic state variables based on specified values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) from which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above
        :return: nothing
        For detailed docstrings on fn functions see calculatestate function, they are the same.
        """
        # Step 1: read which properties are being specified from the combo boxes
        SP = [self._cmb_Property1_2.currentText()[-2:-1], self._cmb_Property2_2.currentText()[-2:-1]]
        if SP[0] == SP[1]:
            self._lbl_Warning.setText("Warning:  You cannot specify the same property twice.")
        else:
            self._lbl_Warning.setText("")
        SP[0] = SP[0].lower()
        SP[1] = SP[1].lower()

        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if SP[0] == 'p' or SP[1] == 'p':
            oFlipped = SP[0] != 'p'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 1:  pt or tp
            if SP1 == 't':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2
                self.t = f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.tSat = round(self.tSat, 3)  # I will compare at 3 three decimal places
                # compare T to TSat
                if self.t < self.tSat or self.t > self.tSat:
                    self.region = "sub-cooled liquid" if self.t < self.tSat else "super-heated vapor"
                    self.makeLabel_1Phase2()
                else:  # this is ambiguous since at saturated temperature
                    self.region = "two-phase"
                    self.stProps = "Region = {:}".format(self.region)
                    self.stProps += "\nPressure = {:0.3f} ({:})".format(self.p, self.p_Units)
                    self.stProps += "\nTemperature = {:0.3f} ({:}) [TSat={:0.3f} ({:})]".format(self.t, self.t_Units,
                                                                                                self.tSat,
                                                                                                self.t_Units)
                    self.stProps += "\nInternal Energy = uf={:0.3f}, ug={:0.3f} ({:})".format(self.uf, self.ug,
                                                                                              self.u_Units)
                    self.stProps += "\nEnthalpy = hf={:0.3f}, hg={:0.3f} ({:})".format(self.hf, self.hg, self.h_Units)
                    self.stProps += "\nEntropy = sf={:0.3f}, sg={:0.3f} ({:})".format(self.sf, self.sg, self.s_Units)
                    self.stProps += "\nSpecific Volume = vf={:0.3f}, vg={:0.3f} ({:})".format(self.vf, self.vg,
                                                                                              self.v_Units)
                    self.stProps += "\nQuality = unknown"
            # case 2: pv or vp
            elif SP1 == 'v':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2
                self.v = f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find T
                    dt = 1.0 if self.v > self.vg else -1.0
                    fn1 = lambda T: self.v - self.steamTable.v_pt(self.p, T)
                    self.t = fsolve(fn1, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.t = self.tSat
                    self.makeLabel_2Phase2()
            # case 3 pu or up
            elif SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2;
                self.u = f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find T
                    dt = 1.0 if self.u > self.ug else -1.0
                    fn3 = lambda T: self.u - self.steamTable.u_pt(self.p, T)
                    self.t = fsolve(fn3, [self.tSat + dt])[0]
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.t = self.tSat
                    self.makeLabel_2Phase2()
            # case 4 ph or hp
            elif SP1 == 'h':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2;
                self.h = f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.t = self.steamTable.t_ph(self.p, self.h)
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.t = self.tSat
                    self.makeLabel_2Phase2()
            # case 5 ps or sp
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2;
                self.s = f2 if not oFlipped else f1;
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.t = self.steamTable.t_ps(self.p, self.s)
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase2()
            # case 6 px or xp
            elif SP1 == 'x':
                self.region = "two-phase"
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.p = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.getSatProps_p(self.p)
                self.t = self.tSat
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.makeLabel_2Phase2()
        elif SP[0] == 't' or SP[1] == 't':
            oFlipped = SP[0] != 't'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 7:  tv or vt
            if SP1 == 'v':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.v = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find P
                    dp = -0.1 if self.v > self.vg else 0.1
                    fn3 = lambda P: [self.v - self.steamTable.v_pt(P, self.t)]
                    self.p = fsolve(fn3, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.p = self.pSat
                    self.makeLabel_2Phase2()
            # case 8:  tu or ut
            elif SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find P
                    dp = 0.1 if self.u > self.ug else -0.1
                    fn8 = lambda P: self.u - self.steamTable.u_pt(self.t, P)
                    self.p = fsolve(fn8, [self.pSat + dp])[0]
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.makeLabel_2Phase2()
            # case 9:  th or ht
            elif SP1 == 'h':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare h to hf and hg
                if self.h < self.hf or self.h > self.hg:
                    self.region = "sub-cooled liquid" if self.h < self.hf else "super-heated vapor"
                    self.p = self.steamTable.p_th(self.t, self.h)
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.x = (self.h - self.hf) / (self.hgf)
                    self.makeLabel_2Phase2()
            # case 10:  ts or st
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1
                self.getSatProps_t(self.t)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    self.p = self.steamTable.p_ts(self.t, self.s)
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase2()
            # case 11:  tx or xt
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.t = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"
                self.getSatProps_t(self.t)
                self.p = self.pSat
                self.x = float(self._le_Property2_2.text())
                self.x = self.clamp(self.x, 0.0, 1.0)
                self.makeLabel_2Phase2()
        elif SP[0] == 'v' or SP[1] == 'v':
            oFlipped = SP[0] != 'v'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 12:  vh or hv
            if SP1 == 'h':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.h = f2 if not oFlipped else f1

                def fn12(P):
                    # could be single phase or two-phase, but both v&h have to match at same x
                    self.getSatProps_p(P)
                    if self.between(self.h, self.hf, self.hg):
                        self.x = (self.h - self.hf) / self.hgf
                        return self.v - (self.vf + self.x * self.vgf)
                    # could be single phase
                    return self.v - self.steamTable.v_ph(P, self.h)

                self.p = fsolve(fn12, [1.0])[0]
                self.t = self.steamTable.t_ph(self.p, self.h)
                self.getSatProps_p(self.p)
                self.vf = round(self.vf, 5)
                self.vg = round(self.vg, 3)
                # compare v to vf and vg
                if self.v < self.vf or self.v > self.vg:
                    self.region = "sub-cooled liquid" if self.v < self.vf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.v - self.vf) / (self.vgf)
                    self.makeLabel_2Phase2()
            # case 13:  vu or uv
            elif SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1

                # use fsolve to fing P&T at this v & u
                def fn13(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u, self.uf, self.ug):
                        self.t = self.tSat
                        self.x = (self.u - self.uf) / self.ugf
                        return [self.v - (self.vf + self.x * self.vgf), 0]
                    return [self.v - self.steamTable.v_pt(PT[0], PT[1]), self.u - self.steamTable.u_pt(PT[0], PT[1])]

                props = fsolve(fn13, [1, 100])
                self.p = props[0]
                self.t = props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t = self.tSat
                    self.makeLabel_2Phase2()
            # case 14:  vs os sv
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1

                def fn14(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.v - self.vf - self.x * self.vgf, 0.0]
                    return [self.v - self.steamTable.v_pt(PT[0], PT[1]),
                            self.s - self.steamTable.s_pt(PT[0], PT[1])]

                props = fsolve(fn14, [1, 100])
                self.p = props[0]
                self.t = props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase2()
            # case 15:  vx or xv
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"

                def fn15(p):
                    self.getSatProps_p(p)
                    return self.v - (self.vf + self.x * self.vgf)

                props = fsolve(fn15, [1])
                self.p = props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()), 0.0, 1.0)
                self.makeLabel_2Phase2()
        elif SP[0] == 'h' or SP[1] == 'h':
            oFlipped = SP[0] != 'h'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 16:  hu or uh
            if SP1 == 'u':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.h = f1 if not oFlipped else f2
                self.u = f2 if not oFlipped else f1

                # use fsolve to fing P&T at this v & u
                def fn16(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.u, self.uf, self.ug):
                        self.x = (self.u - self.uf) / self.ugf
                        return [self.h - self.hf - self.x * self.hgf, 0.0]
                    return [self.h - self.steamTable.h_pt(PT[0], PT[1]), self.u - self.steamTable.u_pt(PT[0], PT[1])]

                props = fsolve(fn16, [1, 100])
                self.p = props[0]
                self.t = props[1]
                self.getSatProps_p(self.p)
                # compare u to uf and ug
                if self.u < self.uf or self.u > self.ug:
                    self.region = "sub-cooled liquid" if self.u < self.uf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.x = (self.u - self.uf) / (self.ugf)
                    self.p = self.pSat
                    self.t = self.tSat
                    self.makeLabel_2Phase2()
            # case 17:  hs or sh
            elif SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.h = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1

                def fn17(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.h - self.hf - self.x * self.hgf, 0.0]
                    return [self.h - self.steamTable.h_pt(PT[0], PT[1]), self.s - self.steamTable.s_pt(PT[0], PT[1])]

                props = fsolve(fn17, [1, 100])
                self.p = props[0]
                self.t = props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase2()
            # case 18:  hx or xh
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.v = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"

                def fn18(p):
                    self.getSatProps_p(p)
                    return self.h - (self.hf + self.x * self.hgf)

                props = fsolve(fn18, [1])
                self.p = props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()), 0.0, 1.0)
                self.makeLabel_2Phase2()
        elif SP[0] == 'u' or SP[1] == 'u':
            oFlipped = SP[0] != 'u'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 19:  us or su
            if SP1 == 's':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.u = f1 if not oFlipped else f2
                self.s = f2 if not oFlipped else f1

                def fn19(PT):
                    self.getSatProps_p(PT[0])
                    if self.between(self.s, self.sf, self.sg):
                        self.x = (self.s - self.sf) / self.sgf
                        return [self.u - self.uf - self.x * self.ugf, 0.0]
                    return [self.u - self.steamTable.u_pt(PT[0], PT[1]),
                            self.s - self.steamTable.s_pt(PT[0], PT[1])]

                props = fsolve(fn19, [1, 100])
                self.p = props[0]
                self.t = props[1]
                self.getSatProps_p(self.p)
                # compare s to sf and sg
                if self.s < self.sf or self.s > self.sg:
                    self.region = "sub-cooled liquid" if self.s < self.sf else "super-heated vapor"
                    # now use P and T
                    self.makeLabel_1Phase2()
                else:  # two-phase
                    self.region = "two-phase"
                    # first calculate quality
                    self.p = self.pSat
                    self.t = self.tSat
                    self.x = (self.s - self.sf) / (self.sgf)
                    self.makeLabel_2Phase2()
            # case 20:  ux or xu
            elif SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.u = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"

                def fn20(p):
                    self.getSatProps_p(p)
                    return self.h - (self.hf + self.x * self.hgf)

                props = fsolve(fn20, [1])
                self.p = props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()), 0.0, 1.0)
                self.makeLabel_2Phase2()
        elif SP[0] == 's' or SP[1] == 's':
            oFlipped = SP[0] != 's'
            SP1 = SP[0] if oFlipped else SP[1]
            # case 21:  sx or xs
            if SP1 == 'x':
                f1 = float(self._le_Property1_2.text())
                f2 = float(self._le_Property2_2.text())
                self.s = f1 if not oFlipped else f2
                self.x = f2 if not oFlipped else f1
                self.region = "two-phase"

                def fn21(p):
                    self.getSatProps_p(p)
                    return self.h - (self.hf + self.x * self.hgf)

                props = fsolve(fn21, [1])
                self.p = props[0]
                self.p = self.pSat
                self.t = self.tSat
                self.x = self.clamp(float(self._le_Property2_2.text()), 0.0, 1.0)
                self.makeLabel_2Phase2()
        self._lbl_StateProperties_2.setText(self.stProps)


#endregion

#region function definitions
def main():
    """
    Entry point function to run the application.

    This function initializes the Qt application if it doesn't already exist,
    connects the aboutToQuit signal to deleteLater to ensure proper cleanup,
    creates the main window, and starts the application event loop.
    """
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())
    pass
#end region

#region function calls
if __name__=="__main__":
    main()
#endregion
