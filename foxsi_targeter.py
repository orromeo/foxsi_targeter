# -*- coding: utf-8 -*-
"""
AUTHOR:  Orlando Romeo, Dan Ryan
UPDATE:  04/06/2024
PURPOSE:
    Find targeting of possible AR locations for FOXSI-4 Campaign
"""
##############################################################################
###############################################################################
# --------------------------------- PyQT Modules
from datetime import datetime
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QTextEdit, \
                            QPushButton, QTableView, QHBoxLayout,QComboBox, \
                            QMessageBox,QDateTimeEdit, QFileDialog,QHeaderView
from PyQt5.QtGui import QStandardItemModel, QStandardItem,QFont
from PyQt5.QtCore import Qt
# --------------------------------- Sunpy/Astropy Modules
import astropy.coordinates
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Helioprojective, RotatedSunFrame
import sunpy.time
import sunpy.map
# --------------------------------- Misc Python Modules
import subprocess
import sys
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
###############################################################################
###############################################################################
# Class for Main GUI
class FOXSITargetGUI(QWidget):
    ###########################################################################
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FOXSI Targeting") # Main GUI
        self.new_window = None                      # Initial SPARCS GUI
        self.initUI()
    ###########################################################################
    def initUI(self):
        # Set Layout
        layout = QVBoxLayout()
        # Horizontal layout for text edit and date input
        hbox = QHBoxLayout()
        #----------------------------------------------------------------------
        # Set title
        title_label = QLabel('FOXSI Active Region Targeting')
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 56px; font-weight: bold;")
        layout.addWidget(title_label)
        # List input
        data_label = QLabel("AR List:")
        data_label.setStyleSheet("font-weight: bold;")
        hbox.addWidget(data_label)
        self.data_text = QTextEdit()
        self.data_text.setFixedHeight(400)
        newline='\n------------------------------------------------------\n'
        self.data_text.setText("Targets Observed: 2024/04/06T00:00:00"+newline) # Set default text
        hbox.addWidget(self.data_text)
        # Launch Date input
        date_label = QLabel("Launch Start Window:\n(LWT)")
        date_label.setAlignment(Qt.AlignCenter)
        date_label.setStyleSheet("font-weight: bold;")
        hbox.addWidget(date_label)        
        ldate = datetime.now()
        self.date_time_edit = QDateTimeEdit(ldate.replace(hour=19,minute=53,second=0,microsecond=0))
        self.date_time_edit.setDisplayFormat("yyyy/MM/ddTHH:mm:ss")  # Set display format
        hbox.addWidget(self.date_time_edit)
        #----------------------------------------------------------------------
        # Add Hbox to Layout
        layout.addLayout(hbox)
        # Create layout for ar locations, and hmi image
        button_layout = QHBoxLayout()
        # Add Science AR Locations Button
        convert_button = QPushButton("SCIENCE: Targeting")
        convert_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert_button.clicked.connect(self.science_target)
        button_layout.addWidget(convert_button)
        # Add Science AR Locations Button
        convert2_button = QPushButton("SPARCS: Targeting")
        convert2_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert2_button.clicked.connect(self.sparcs_target)
        button_layout.addWidget(convert2_button)
        # Add HMI Button
        picture_button = QPushButton("SAAS: HMI Image")
        picture_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        picture_button.clicked.connect(self.show_hmi)
        button_layout.addWidget(picture_button)
        # Add buttons to layout
        layout.addLayout(button_layout)
        #----------------------------------------------------------------------
        # Science Table Output
        space_label = QLabel(" ")
        space_label.setAlignment(Qt.AlignCenter)
        space_label.setFixedHeight(100)
        space_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addWidget(space_label)
        table_label = QLabel("FOXSI SCIENCE: AR Targeting Table")
        table_label.setAlignment(Qt.AlignCenter)
        table_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addWidget(table_label)
        self.table_view = QTableView()
        layout.addWidget(self.table_view)
        #----------------------------------------------------------------------
        # Export to CSV
        convert_button = QPushButton("Export to CSV")
        convert_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert_button.clicked.connect(lambda: self.science_exportToCSV(self.date_time_edit.dateTime().toString("yyyy_MM_dd")))
        layout.addWidget(convert_button)
        #----------------------------------------------------------------------
        # Set Layout
        self.setLayout(layout)
    ###########################################################################
    # Function that contracts off disk AR to right on the limb
    def contract_to_disk(self,Tx,Ty,factor=1):
        rsun_arcsec = 959.63
        Tx_unit     = Tx/np.sqrt(Tx**2+Ty**2)
        Ty_unit     = Ty/np.sqrt(Tx**2+Ty**2)
        Tx_onlimb   = Tx_unit*rsun_arcsec*factor # Factor of solar radii to consider
        Ty_onlimb   = Ty_unit*rsun_arcsec*factor
        return Tx_onlimb,Ty_onlimb
    ###########################################################################
    # Transform AR Location to FOXSI Observer during Launch Window
    def transform_target_to_foxsi_hpc(self,Tx,Ty, obstime,launch_starttime, **kwargs):
        # Get launch window times        
        foxsi_times = sunpy.time.parse_time(launch_starttime) + np.arange(0, 5, 1) * u.hour
        # Get observed time of AR
        obs_time   = sunpy.time.parse_time(obstime,scale='utc')
        #----------------------------------------------------------------------
        # Get coordinate of FOXSI at apogee.
        poker_lon, poker_lat, poker_altitude = -147.47*u.deg, 65.12*u.deg, 197*u.m
        apogee = 400*u.km
        foxsi_loc = astropy.coordinates.EarthLocation(lon=poker_lon, lat=poker_lat, height=poker_altitude + apogee)
        foxsi_skycoord = SkyCoord(foxsi_loc.get_gcrs(foxsi_times))
        # Find observer of AR
        target_observer = astropy.coordinates.get_body("Earth", time=obs_time)
        #----------------------------------------------------------------------
        # Iterate until AR is on disk in case it is off disk
        new_Tx = np.nan
        factor = 1
        while np.any(np.isnan(new_Tx)):
            # Create Helioprojective coordinates
            planning_target = SkyCoord(Tx * u.arcsec, Ty * u.arcsec,
                              frame=Helioprojective(observer=target_observer))
            # Differentially rotate target to launch time.
            target_diffrot = SkyCoord(\
                             RotatedSunFrame(base=planning_target.frame,\
                                             rotated_time=foxsi_times,\
                                             rotation_model=kwargs.get("rotation_model", "howard")))
            # Transform to FOXSI observer view.
            flight_target = target_diffrot.transform_to(Helioprojective(observer=foxsi_skycoord))
            new_Tx = flight_target.Tx
            # Check if nan due to off limb AR
            if any(np.isnan(new_Tx)):
                Tx,Ty   = self.contract_to_disk(Tx,Ty,factor=factor)
                factor -= 0.001
        #----------------------------------------------------------------------
        return Tx,Ty,flight_target.Tx,flight_target.Ty
    ###########################################################################
    def science_target(self):
        #----------------------------------------------------------------------
        # Get the input data
        data = self.data_text.toPlainText().strip()
        # Check if data is empty
        if not data:
            QMessageBox.warning(self, "Error", "Please Enter AR Data.")
            return
        # Split the data into rows
        rows = data.split('\n')
        # Check for AR Data
        if len(rows) < 3:
            QMessageBox.warning(self, "Error", "Please Enter AR Data.")
            return
        #----------------------------------------------------------------------
        # Get AR data
        rows_data = rows[2:]
        # Get AR Observation Date
        obsdate = (rows[0].split(' '))[-1]
        # Get the selected launch date
        launchdate = self.date_time_edit.dateTime().toString("yyyy/MM/ddTHH:mm:ss")
        # Find time difference between observation and launch
        rot_duration = (datetime.strptime(obsdate, "%Y/%m/%dT%H:%M:%S") - \
                        datetime.strptime(launchdate, "%Y/%m/%dT%H:%M:%S")).total_seconds()/3600
        # Initialize Table contents
        table_data   = []
        row_names    = []
        column_names = []
        #----------------------------------------------------------------------
        # Iterate each AR
        for ri,row in enumerate(rows_data):
            parts = row.split(': ')
            if len(parts) < 2:
                QMessageBox.warning(self, "Error", "Invalid data format.")
                return
            row_names.append(" "+str(ri+1)+". "+parts[0]+"         ")
            values = parts[1].replace(' (', ' ').replace(')', '').replace(',', ' ').split()
            fvals  = []
            # Assign HPC Values
            Tx = float(values[-2].replace('"', ''))
            Ty = float(values[-1].replace('"', ''))
            newTx,newTy,launchlon,launchlat = self.transform_target_to_foxsi_hpc(Tx,Ty,obsdate,launchdate)
            # Check if on solar disk
            if newTx == Tx:
                fvals.append('Yes')
            else:
                fvals.append('No')
            fvals.append("({:+.2f}, {:+.2f})".format(newTx, newTy))
            #------------------------------------------------------------------
            launchcoor = ["({:+.2f}, {:+.2f})".format(x, y) \
                          for x, y, in zip(launchlon.to_value(u.arcsec),launchlat.to_value(u.arcsec))]
            fvals.extend(launchcoor)
            # Assign Column Values in Row
            offsetlaunchcoor = ["({:+.2f}, {:+.2f})".format(x+149.21, y-156.15) \
                          for x,y in zip(launchlon.to_value(u.arcsec),launchlat.to_value(u.arcsec))]
            # Add Column Values
            fvals.extend(offsetlaunchcoor)
            table_data.append(fvals)
        #----------------------------------------------------------------------
        # Assign Column Names
        column_names = ["On Disk?",' \nSDO\n(θx,θy) ["]\nLWT{:+.1f}hr\n   '.format(rot_duration)]
        launchnames  =  [f' \nOPTICS\n(θx,θy) ["]\nLWT+{i}hr\n   ' for i in range(0, len(launchlon))]
        offsetnames  =  [f' \nSPARCS\n(θx,θy) ["]\nLWT+{i}hr\n   ' for i in range(0, len(launchlon))]
        column_names.extend(launchnames)
        column_names.extend(offsetnames)
        # Check for Errors
        if not row_names or not column_names or not table_data:
            QMessageBox.warning(self, "Error", "Data could not be processed.")
            return
        #----------------------------------------------------------------------
        # Set up table
        model = QStandardItemModel(len(row_names), len(column_names))
        model.setHorizontalHeaderLabels(column_names)
        model.setVerticalHeaderLabels(row_names)
        for row_idx, row in enumerate(table_data):
            # Extract every other value as it corresponds to the data we need
            for col_idx, value in enumerate(row):
                item = QStandardItem(value)
                item.setTextAlignment(Qt.AlignCenter)
                model.setItem(row_idx, col_idx, item)
        # Set font style for row names and column names
        font = QFont()
        font.setBold(True)
        for row_idx in range(len(row_names)):
            model.setHeaderData(row_idx, Qt.Vertical, font, role=Qt.FontRole)
        for col_idx in range(len(column_names)):
            model.setHeaderData(col_idx, Qt.Horizontal, font, role=Qt.FontRole)
        #----------------------------------------------------------------------
        # Set gridstyle
        self.table_view.setStyleSheet(
            "QTableView {gridline-color: black;}"
            "QHeaderView::section {background-color: white; border: 2px solid black;}")
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # Stretch Vertical Cells
        vheader = self.table_view.verticalHeader()  # Corrected attribute name
        vheader.setSectionResizeMode(QHeaderView.Stretch)  # Stretch columns to fill available space    
        # Show Grid
        self.table_view.setShowGrid(True)
        # Display the model in the table view
        self.table_view.setModel(model)
        
    ###########################################################################
    def sparcs_target(self):
        #----------------------------------------------------------------------
        # Get the input data
        data = self.data_text.toPlainText().strip()
        # Check if data is empty
        if not data:
            QMessageBox.warning(self, "Error", "Please Enter AR Data.")
            return
        # Split the data into rows
        rows = data.split('\n')
        # Check for AR Data
        if len(rows) < 3:
            QMessageBox.warning(self, "Error", "Please Enter AR Data.")
            return
        #----------------------------------------------------------------------
        # Create new GUI
        self.new_window = SPARCSGUI(self)
        self.new_window.showMaximized()
    ###########################################################################
    def show_hmi(self):
        # Open a file explorer dialog to select a file
        default_dir = os.path.join(os.getcwd(), "HMI")
        img_dir = os.path.join(default_dir, "Images")
        data_dir = os.path.join(default_dir, "Data")
        # Check if the directory exists, if not create it
        if not os.path.exists(img_dir):
            os.makedirs(img_dir)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Picture",data_dir, "Image Files (*.png *.jpg *.jpeg *.bmp *.fits)")
        # Check if a file is selected
        if file_path:
            hmi_map = sunpy.map.Map(file_path)
            hmi_rotated = hmi_map.rotate(angle=-90*u.degree,order=3)
            # aia_rotated = hmi_map.rotate(angle=-90*u.degree,order=3)
            # fig = plt.figure()
            # ax = fig.add_subplot(projection=aia_rotated)
            # aia_rotated.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
            # aia_rotated.draw_limb(axes=ax)
            # aia_rotated.draw_grid(axes=ax)
            # plt.show()
            #hmi_rotated.plot(axes=ax)
            # Plot Data
            data = ((hmi_rotated.data))
            data[np.where(np.isnan(data))] = 0
            datasz = np.shape(data)
            fig = plt.figure(figsize=(20,20))
            ax = fig.add_subplot(projection=hmi_rotated)
            ax.imshow(data, cmap='hot',extent=[0,datasz[0],0,datasz[1]])
            # Create a custom ticker formatter
            if 0:
                scl = hmi_rotated.scale[0].value
                label_format = '{:,.0f}'
                ticks_loc = ax.get_xticks().tolist()
                ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
                ax.set_xticklabels([label_format.format(x) for x in ((ax.get_xticks()-datasz[0]/2)*scl)])
                ticks_loc = ax.get_yticks().tolist()
                ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
                ax.set_yticklabels([label_format.format(x) for x in ((ax.get_yticks()-datasz[0]/2)*scl)])
            ax.set_title("SAAS: HMI IMAGE", fontsize=38)
            plt.show()
            # Save and Open File
            sfile = os.path.join(img_dir, "HMI_IMAGE_"+hmi_map.date.strftime('%Y-%m-%dT%H%M%S'))
            plt.savefig(sfile)
            plt.close()
            subprocess.Popen(sfile+".png", shell=True)
        else:
            # Display the selected picture (placeholder code)
            QMessageBox.warning(self, "Error", "No Image!")
    ###########################################################################
    def science_exportToCSV(self,launchdate):
        # Create filename
        filedir  = os.path.join(os.getcwd(), "Science-Target_Tables")
        # Check if the directory exists, if not create it
        if not os.path.exists(filedir):
            os.makedirs(filedir)
        filename = os.path.join(filedir, "Science-Target_Tables_"+launchdate+".csv")
        # Access the table model
        model = self.table_view.model()
        if model is not None:
            # Open the CSV file for writing
            with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                # Write column headers
                headers = ['Target']
                headers.extend([model.headerData(col, Qt.Horizontal).replace("OPTICS","FOXSI "
                                                                    ).replace("\n",""
                                                                    ).replace("Δ"," del"
                                                                    ).replace("θ"," T"
                                                                    ).replace("θ"," T"
                                                                    ).replace("hr","hr)"
                                                                    ).replace("LWT"," (LWT") for col in range(model.columnCount()-5)])
                writer.writerow(headers)
                # Write data rows
                for row in range(model.rowCount()):
                    row_name = model.headerData(row, Qt.Vertical)
                    row_data = [row_name] + [model.index(row, col).data() for col in range(model.columnCount()-5)]
                    writer.writerow(row_data)    
            QMessageBox.information(self, 'CSV Generated', 'CSV File Created!')
        else:
            QMessageBox.warning(self, "Error", "No Data!")
            return
    ###########################################################################
    def closeEvent(self,event):
        QApplication.quit()      
###############################################################################
###############################################################################            
class SPARCSGUI(QWidget):
    ###########################################################################
    def __init__(self, parent):
        super().__init__()
        self.setWindowTitle("SPARCS Targeting")
        self.initUI(parent)
    ###########################################################################
    def initUI(self,parent):
        # Create Layout
        layout = QVBoxLayout()
        # Horizontal layout for text edit and date input
        hbox = QHBoxLayout()
        title_label = QLabel("SPARCS Targeting")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addLayout(hbox)
        # Add title label
        title_label = QLabel("SPARCS Targeting")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addWidget(title_label)
        # Create layout for "Aligned" label and combo box
        aligned_layout = QHBoxLayout()
        # Add label for "Aligned" text
        aligned_label = QLabel("Solar-North Alignment:")
        aligned_label.setAlignment(Qt.AlignCenter)
        aligned_label.setFixedHeight(200)
        aligned_label.setStyleSheet("font-size: 38px; font-weight: bold;")
        aligned_layout.addWidget(aligned_label)
        # Create combo box for drop-down selection
        self.comboBox = QComboBox()
        self.comboBox.addItems(["Payload +0", "Payload +180"])
        self.comboBox.setStyleSheet("font-size: 38px; font-weight: bold;")
        aligned_layout.addWidget(self.comboBox)
        self.comboBox.currentIndexChanged.connect(self.multiplyRowValues)
        # Add aligned layout to main layout
        layout.addLayout(aligned_layout)
        #----------------------------------------------------------------------
        # Initialize Table
        self.table_view = QTableView()
        layout.addWidget(self.table_view)
        self.setLayout(layout)
        #----------------------------------------------------------------------
        # Get the input data
        data = parent.data_text.toPlainText().strip()
        # Split the data into rows
        rows = data.split('\n')
        # Get AR data
        rows_data = rows[2:]
        # Get AR Observation Date
        obsdate = (rows[0].split(' '))[-1]
        # Get the selected launch date
        launchdate = parent.date_time_edit.dateTime().toString("yyyy/MM/ddTHH:mm:ss")
        # Initialize Table contents
        table_data   = []
        row_names    = []
        column_names = []
        #----------------------------------------------------------------------
        # Iterate each AR
        for ri,row in enumerate(rows_data):
            parts = row.split(': ')
            if len(parts) < 2:
                QMessageBox.warning(self, "Error", "Invalid data format.")
                return
            row_names.append(" "+str(ri+1)+". "+parts[0]+"         ")
            values = parts[1].replace(' (', ' ').replace(')', '').replace(',', ' ').split()
            fvals  = []
            # Assign HPC Values
            Tx = float(values[-2].replace('"', ''))
            Ty = float(values[-1].replace('"', ''))
            newTx,newTy,launchlon,launchlat = parent.transform_target_to_foxsi_hpc(Tx,Ty,obsdate,launchdate)
            #------------------------------------------------------------------
            # Check payload alignment
            index = self.comboBox.currentIndex()
            if self.comboBox.itemText(index) == "Payload +0":
                multiplier = 1
            else:
                multiplier= -1
            # Compute new delta values
            lonround = np.round(launchlon.to_value(u.arcsec)+149.21)
            latround = np.round(launchlat.to_value(u.arcsec)-156.15)
            delx = 10
            dely = 90
            lonsprc = [multiplier*np.round((x - lonround[2]) / delx) * delx for x in lonround]
            latsprc = [multiplier*np.round((x - latround[2]) / dely) * dely for x in latround]
            lonsprc[2] = lonround[2]
            latsprc[2] = latround[2]
            # Assign Column Values in Row
            offsetlaunchcoor = ["{:+}".format(x) for x in lonsprc]
            offsetlaunchcoor[2] = "(  {:+.2f}, {:+.2f}  )".format(lonround[2],latround[2])
            # Add Column Values
            fvals.extend(offsetlaunchcoor)
            table_data.append(fvals)
        #----------------------------------------------------------------------
        # Assign Column Names
        column_names = [f' \nSPARCS\nΔθx ["]\nLWT+{i}hr\n   ' for i in range(0,int( (len(launchlon)-1)/2))]
        column_names.extend([' \nSPARCS\n( θx, θy ) ["]\nLWT+2hr\n   '])
        column_names.extend([f' \nSPARCS\nΔθx ["]\nLWT+{i}hr\n   ' for i in range(3,len(launchlon))])
        # Check for Errors
        if not row_names or not column_names or not table_data:
            QMessageBox.warning(self, "Error", "Data could not be processed.")
            return
        #----------------------------------------------------------------------
        # Set up table
        model = QStandardItemModel(len(row_names), len(column_names))
        model.setHorizontalHeaderLabels(column_names)
        model.setVerticalHeaderLabels(row_names)
        for row_idx, row in enumerate(table_data):
            # Extract every other value as it corresponds to the data we need
            for col_idx, value in enumerate(row):
                item = QStandardItem(value)
                item.setTextAlignment(Qt.AlignCenter)
                model.setItem(row_idx, col_idx, item)
        # Set font style for row names and column names
        font = QFont()
        font.setBold(True)
        for row_idx in range(len(row_names)):
            model.setHeaderData(row_idx, Qt.Vertical, font, role=Qt.FontRole)
        for col_idx in range(len(column_names)):
            model.setHeaderData(col_idx, Qt.Horizontal, font, role=Qt.FontRole)
        #----------------------------------------------------------------------
        # Set gridstyle
        self.table_view.setStyleSheet(
            "QTableView {gridline-color: black;}"
            "QHeaderView::section {background-color: white; border: 2px solid black;}")
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        # Stretch Vertical Cells
        vheader = self.table_view.verticalHeader()  # Corrected attribute name
        vheader.setSectionResizeMode(QHeaderView.Stretch)  # Stretch columns to fill available space    
        # Show Grid
        self.table_view.setShowGrid(True)
        # Display the model in the table view
        self.table_view.setModel(model)   
        #----------------------------------------------------------------------
        # Export to CSV
        convert_button = QPushButton("Export to CSV")
        convert_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert_button.clicked.connect(lambda: self.sparcs_exportToCSV(parent.date_time_edit.dateTime().toString("yyyy_MM_dd")))
        layout.addWidget(convert_button)
    ###########################################################################
    def multiplyRowValues(self, index):
        # Access the table model
        model = self.table_view.model()
        # Iterate through each cell in the table
        for row in range(model.rowCount()):
            for col in range(model.columnCount()):
                # Get the current cell value
                cell_value = model.index(row, col).data()
                # Convert the cell value to a float (assuming the value is convertible to a number)
                try:
                    cell_value = float(cell_value)
                except ValueError:
                    continue  # Skip non-numeric values
                # Multiply the cell value by -1 based on the selected index
                #if index == 0:  # If "Payload +180" is selected
                cell_value *= -1
                # Set the updated value back to the cell
                model.setData(model.index(row, col), str(cell_value))
    ###########################################################################
    def sparcs_exportToCSV(self,launchdate):
        # Create filename
        filedir  = os.path.join(os.getcwd(), "SPARCS-Target_Tables")
        # Check if the directory exists, if not create it
        if not os.path.exists(filedir):
            os.makedirs(filedir)
        filename = os.path.join(filedir, "SPARCS-Target_Tables_"+launchdate+".csv")
        # Access the table model
        model = self.table_view.model()
        # Open the CSV file for writing
        with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            # Write column headers
            headers = ['Target']
            headers.extend([model.headerData(col, Qt.Horizontal).replace("SPARCS",""
                                                               ).replace("\n",""
                                                               ).replace("Δ"," del"
                                                               ).replace("θ"," T"
                                                               ).replace("θ"," T"
                                                               ).replace("hr","hr)"
                                                               ).replace("LWT"," (LWT") for col in range(model.columnCount())])
            writer.writerow(headers)
            # Write data rows
            for row in range(model.rowCount()):
                row_name = model.headerData(row, Qt.Vertical)
                row_data = [row_name] + [model.index(row, col).data() for col in range(model.columnCount())]
                writer.writerow(row_data)
        QMessageBox.information(self, 'CSV Generated', 'CSV File Created!')
###############################################################################
###############################################################################
# Main 
if __name__ == '__main__':
    app       = QApplication(sys.argv)
    converter = FOXSITargetGUI()
    # Maximize the window
    converter.showMaximized()  
    sys.exit(app.exec_())