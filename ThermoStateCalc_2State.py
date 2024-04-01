# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ThermoStateCalc5_2State.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui__frm_StateCalculator(object):
    def setupUi(self, _frm_StateCalculator):
        _frm_StateCalculator.setObjectName("_frm_StateCalculator")
        _frm_StateCalculator.resize(702, 459)
        self.verticalLayout = QtWidgets.QVBoxLayout(_frm_StateCalculator)
        self.verticalLayout.setObjectName("verticalLayout")
        self._grp_Units = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_Units.setObjectName("_grp_Units")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self._grp_Units)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self._rdo_SI = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_SI.setChecked(True)
        self._rdo_SI.setObjectName("_rdo_SI")
        self.horizontalLayout_2.addWidget(self._rdo_SI)
        self._rdo_English = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_English.setObjectName("_rdo_English")
        self.horizontalLayout_2.addWidget(self._rdo_English)
        self.verticalLayout.addWidget(self._grp_Units)
        self._grp_SpecifiedProperties = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_SpecifiedProperties.setObjectName("_grp_SpecifiedProperties")
        self.gridLayout = QtWidgets.QGridLayout(self._grp_SpecifiedProperties)
        self.gridLayout.setObjectName("gridLayout")
        self._pb_Calculate = QtWidgets.QPushButton(self._grp_SpecifiedProperties)
        self._pb_Calculate.setObjectName("_pb_Calculate")
        self.gridLayout.addWidget(self._pb_Calculate, 5, 0, 1, 1)
        self._le_Property2_2 = QtWidgets.QLineEdit(self._grp_SpecifiedProperties)
        self._le_Property2_2.setObjectName("_le_Property2_2")
        self.gridLayout.addWidget(self._le_Property2_2, 3, 9, 1, 1)
        self._lbl_Property1_2 = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property1_2.setObjectName("_lbl_Property1_2")
        self.gridLayout.addWidget(self._lbl_Property1_2, 1, 6, 1, 1)
        self._le_Property2 = QtWidgets.QLineEdit(self._grp_SpecifiedProperties)
        self._le_Property2.setObjectName("_le_Property2")
        self.gridLayout.addWidget(self._le_Property2, 3, 3, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 3, 5, 1, 1)
        self._cmb_Property2 = QtWidgets.QComboBox(self._grp_SpecifiedProperties)
        self._cmb_Property2.setObjectName("_cmb_Property2")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self.gridLayout.addWidget(self._cmb_Property2, 2, 3, 1, 2)
        self.label_2 = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 6, 1, 1)
        self._le_Property1_2 = QtWidgets.QLineEdit(self._grp_SpecifiedProperties)
        self._le_Property1_2.setObjectName("_le_Property1_2")
        self.gridLayout.addWidget(self._le_Property1_2, 3, 6, 1, 1)
        self._le_Property1 = QtWidgets.QLineEdit(self._grp_SpecifiedProperties)
        self._le_Property1.setObjectName("_le_Property1")
        self.gridLayout.addWidget(self._le_Property1, 3, 0, 1, 1)
        self._lbl_Property1 = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property1.setObjectName("_lbl_Property1")
        self.gridLayout.addWidget(self._lbl_Property1, 1, 0, 1, 1)
        self._cmb_Property2_2 = QtWidgets.QComboBox(self._grp_SpecifiedProperties)
        self._cmb_Property2_2.setObjectName("_cmb_Property2_2")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self._cmb_Property2_2.addItem("")
        self.gridLayout.addWidget(self._cmb_Property2_2, 2, 9, 1, 2)
        self._cmb_Property1_2 = QtWidgets.QComboBox(self._grp_SpecifiedProperties)
        self._cmb_Property1_2.setObjectName("_cmb_Property1_2")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self.gridLayout.addWidget(self._cmb_Property1_2, 2, 6, 1, 2)
        spacerItem1 = QtWidgets.QSpacerItem(1, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 3, 8, 1, 1)
        self._lbl_Property1_2_Units = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property1_2_Units.setObjectName("_lbl_Property1_2_Units")
        self.gridLayout.addWidget(self._lbl_Property1_2_Units, 3, 7, 1, 1)
        self._lbl_Property2_2 = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property2_2.setObjectName("_lbl_Property2_2")
        self.gridLayout.addWidget(self._lbl_Property2_2, 1, 9, 1, 1)
        self._lbl_Property1_Units = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property1_Units.setObjectName("_lbl_Property1_Units")
        self.gridLayout.addWidget(self._lbl_Property1_Units, 3, 1, 1, 1)
        self._lbl_Property2 = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property2.setObjectName("_lbl_Property2")
        self.gridLayout.addWidget(self._lbl_Property2, 1, 3, 1, 1)
        self.label = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self._cmb_Property1 = QtWidgets.QComboBox(self._grp_SpecifiedProperties)
        self._cmb_Property1.setObjectName("_cmb_Property1")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self._cmb_Property1.addItem("")
        self.gridLayout.addWidget(self._cmb_Property1, 2, 0, 1, 2)
        self._lbl_Warning = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Warning.setText("")
        self._lbl_Warning.setObjectName("_lbl_Warning")
        self.gridLayout.addWidget(self._lbl_Warning, 4, 0, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(0, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 3, 2, 1, 1)
        self._pb_Calculate_2 = QtWidgets.QPushButton(self._grp_SpecifiedProperties)
        self._pb_Calculate_2.setObjectName("_pb_Calculate_2")
        self.gridLayout.addWidget(self._pb_Calculate_2, 5, 6, 1, 1)
        self._lbl_Property2_Units = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property2_Units.setObjectName("_lbl_Property2_Units")
        self.gridLayout.addWidget(self._lbl_Property2_Units, 3, 4, 1, 1)
        self._lbl_Property2_2_Units = QtWidgets.QLabel(self._grp_SpecifiedProperties)
        self._lbl_Property2_2_Units.setObjectName("_lbl_Property2_2_Units")
        self.gridLayout.addWidget(self._lbl_Property2_2_Units, 3, 10, 1, 1)
        self.verticalLayout.addWidget(self._grp_SpecifiedProperties)
        self._grp_StateProperties = QtWidgets.QGroupBox(_frm_StateCalculator)
        self._grp_StateProperties.setObjectName("_grp_StateProperties")
        self.gridLayout_2 = QtWidgets.QGridLayout(self._grp_StateProperties)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self._lbl_State = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_State.setObjectName("_lbl_State")
        self.gridLayout_2.addWidget(self._lbl_State, 0, 1, 1, 1)
        self._lbl_State_2 = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_State_2.setObjectName("_lbl_State_2")
        self.gridLayout_2.addWidget(self._lbl_State_2, 0, 3, 1, 1)
        self._lbl_StateChange = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_StateChange.setObjectName("_lbl_StateChange")
        self.gridLayout_2.addWidget(self._lbl_StateChange, 0, 5, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem3, 1, 0, 1, 1)
        self._lbl_StateProperties = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_StateProperties.setObjectName("_lbl_StateProperties")
        self.gridLayout_2.addWidget(self._lbl_StateProperties, 1, 1, 1, 1)
        spacerItem4 = QtWidgets.QSpacerItem(65, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem4, 1, 2, 1, 1)
        self._lbl_StateProperties_2 = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_StateProperties_2.setObjectName("_lbl_StateProperties_2")
        self.gridLayout_2.addWidget(self._lbl_StateProperties_2, 1, 3, 1, 1)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem5, 1, 4, 1, 1)
        self._lbl_StateChangeProperties = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_StateChangeProperties.setObjectName("_lbl_StateChangeProperties")
        self.gridLayout_2.addWidget(self._lbl_StateChangeProperties, 1, 5, 1, 1)
        spacerItem6 = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem6, 1, 6, 1, 1)
        self._lbl_SatLiqProps = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_SatLiqProps.setText("")
        self._lbl_SatLiqProps.setObjectName("_lbl_SatLiqProps")
        self.gridLayout_2.addWidget(self._lbl_SatLiqProps, 2, 1, 1, 1)
        self._lbl_SatVapProps = QtWidgets.QLabel(self._grp_StateProperties)
        self._lbl_SatVapProps.setText("")
        self._lbl_SatVapProps.setObjectName("_lbl_SatVapProps")
        self.gridLayout_2.addWidget(self._lbl_SatVapProps, 3, 1, 1, 1)
        self.verticalLayout.addWidget(self._grp_StateProperties)
        spacerItem7 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem7)

        self.retranslateUi(_frm_StateCalculator)
        self._cmb_Property2.setCurrentIndex(1)
        self._cmb_Property2_2.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(_frm_StateCalculator)

    def retranslateUi(self, _frm_StateCalculator):
        _translate = QtCore.QCoreApplication.translate
        _frm_StateCalculator.setWindowTitle(_translate("_frm_StateCalculator", "Thermodynamic State Calculator"))
        self._grp_Units.setTitle(_translate("_frm_StateCalculator", "System of Units"))
        self._rdo_SI.setText(_translate("_frm_StateCalculator", "SI"))
        self._rdo_English.setText(_translate("_frm_StateCalculator", "English"))
        self._grp_SpecifiedProperties.setTitle(_translate("_frm_StateCalculator", "Specified Properties"))
        self._pb_Calculate.setText(_translate("_frm_StateCalculator", "Calculate State 1"))
        self._le_Property2_2.setText(_translate("_frm_StateCalculator", "100.0"))
        self._lbl_Property1_2.setText(_translate("_frm_StateCalculator", "Property 1"))
        self._le_Property2.setText(_translate("_frm_StateCalculator", "100.0"))
        self._cmb_Property2.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property2.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property2.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property2.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property2.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property2.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property2.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self.label_2.setText(_translate("_frm_StateCalculator", "State 2"))
        self._le_Property1_2.setText(_translate("_frm_StateCalculator", "1.0"))
        self._le_Property1.setText(_translate("_frm_StateCalculator", "1.0"))
        self._lbl_Property1.setText(_translate("_frm_StateCalculator", "Property 1"))
        self._cmb_Property2_2.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property2_2.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property2_2.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property2_2.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property2_2.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property2_2.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property2_2.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._cmb_Property1_2.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property1_2.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property1_2.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property1_2.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property1_2.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property1_2.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property1_2.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._lbl_Property1_2_Units.setText(_translate("_frm_StateCalculator", "Bar"))
        self._lbl_Property2_2.setText(_translate("_frm_StateCalculator", "Property 2"))
        self._lbl_Property1_Units.setText(_translate("_frm_StateCalculator", "Bar"))
        self._lbl_Property2.setText(_translate("_frm_StateCalculator", "Property 2"))
        self.label.setText(_translate("_frm_StateCalculator", "State 1"))
        self._cmb_Property1.setItemText(0, _translate("_frm_StateCalculator", "Pressure (p)"))
        self._cmb_Property1.setItemText(1, _translate("_frm_StateCalculator", "Temperature (T)"))
        self._cmb_Property1.setItemText(2, _translate("_frm_StateCalculator", "Quality (x)"))
        self._cmb_Property1.setItemText(3, _translate("_frm_StateCalculator", "Specific Internal Energy (u)"))
        self._cmb_Property1.setItemText(4, _translate("_frm_StateCalculator", "Specific Enthalpy (h)"))
        self._cmb_Property1.setItemText(5, _translate("_frm_StateCalculator", "Specific Volume (v)"))
        self._cmb_Property1.setItemText(6, _translate("_frm_StateCalculator", "Specific Entropy (s)"))
        self._pb_Calculate_2.setText(_translate("_frm_StateCalculator", "Calculate State 2 / State Change"))
        self._lbl_Property2_Units.setText(_translate("_frm_StateCalculator", "C"))
        self._lbl_Property2_2_Units.setText(_translate("_frm_StateCalculator", "C"))
        self._grp_StateProperties.setTitle(_translate("_frm_StateCalculator", "State Properties"))
        self._lbl_State.setText(_translate("_frm_StateCalculator", "State1:  saturated"))
        self._lbl_State_2.setText(_translate("_frm_StateCalculator", "State2:  saturated"))
        self._lbl_StateChange.setText(_translate("_frm_StateCalculator", "State Change:"))
        self._lbl_StateProperties.setText(_translate("_frm_StateCalculator", "Pressure = 1000 kPa\n"
"Temperature = 100 C\n"
"X = 1.0"))
        self._lbl_StateProperties_2.setText(_translate("_frm_StateCalculator", "Pressure = 1000 kPa\n"
"Temperature = 100 C\n"
"X = 1.0"))
        self._lbl_StateChangeProperties.setText(_translate("_frm_StateCalculator", "Pressure = 1000 kPa\n"
"Temperature = 100 C\n"
"X = 1.0"))
