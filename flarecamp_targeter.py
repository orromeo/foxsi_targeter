# -*- coding: utf-8 -*-
"""
AUTHOR:  Orlando Romeo, Dan Ryan
UPDATE:  04/09/2024
PURPOSE:
    Find targeting of possible AR locations for FOXSI-4 and Hi-C Campaign
"""
###############################################################################
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

# Campaign Inputs
###############################################################################
# Set alphabet list
azlist = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
colors = ['red', 'yellow', 'blue', 'green', 'purple','cyan','green','grey','brown','black','silver','white']
# set default launch window time
lwtime     = [19,53,0,0] #Hour, Minute, Second, Milisecond
lwdur      = 4 # Duration of launch window
lcl_time   = -8 # UTC to Local Time
# Create dictionary
foxsi_dict = {'mission':'FOXSI','sparcsoffset':[149.21,-156.15],'dtime':1,'dlon':10,'dlat':90}
hic_dict = {'mission':'Hi-C','sparcsoffset':[np.nan,np.nan],'dtime':9,'dlon':5,'dlat':90}
# Set rocket coordinates
rocket_lon, rocket_lat, rocket_altitude = -147.47*u.deg, 65.12*u.deg, 197*u.m
rocket_apogee = 400*u.km
###############################################################################
###############################################################################
###############################################################################
# Class for Main GUI
class FOXSITargetGUI(QWidget):
    ###########################################################################
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Flare Campaign Targeting") # Main GUI
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
        title_label = QLabel('Flare Campaign: Active Region Targeting')
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
        self.date_time_edit = QDateTimeEdit(ldate.replace(hour=lwtime[0],minute=lwtime[1],second=lwtime[2],microsecond=lwtime[3]))
        self.date_time_edit.setDisplayFormat("yyyy/MM/ddTHH:mm:ss")  # Set display format
        hbox.addWidget(self.date_time_edit)
        # Add Hbox to Layout
        layout.addLayout(hbox)
        #----------------------------------------------------------------------
        # Create layout for ar locations, and hmi image
        button_layout = QHBoxLayout()
        # Add Science AR Locations Button
        convert_button = QPushButton("FOXSI SCIENCE Target Table")
        convert_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert_button.clicked.connect(lambda: self.science_target(foxsi_dict))
        button_layout.addWidget(convert_button)
        # Add Science AR Locations Button
        convert2_button = QPushButton("FOXSI SPARCS Target Table")
        convert2_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert2_button.clicked.connect(lambda: self.sparcs_target(foxsi_dict))
        button_layout.addWidget(convert2_button)
        # Add Science AR Locations Button
        converth_button = QPushButton("Hi-C SCIENCE Target Table")
        converth_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        converth_button.clicked.connect(lambda: self.science_target(hic_dict))
        button_layout.addWidget(converth_button)
        # Add Science AR Locations Button
        convert3_button = QPushButton("Hi-C SPARCS Target Table")
        convert3_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert3_button.clicked.connect(lambda: self.sparcs_target(hic_dict))
        button_layout.addWidget(convert3_button)
        # Add buttons to layout
        layout.addLayout(button_layout)
        #----------------------------------------------------------------------
        # Create layout for ar locations, and hmi image
        im_layout = QHBoxLayout()
        # Add AIA
        aia_button = QPushButton("AIA MAP")
        aia_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        aia_button.clicked.connect(lambda: self.show_map("aia"))
        im_layout.addWidget(aia_button)
        # Add HMI MAG
        hmim_button = QPushButton("HMI MAG MAP")
        hmim_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        hmim_button.clicked.connect(lambda: self.show_map("maghmi"))
        im_layout.addWidget(hmim_button)
        # Add HMI CONT
        hmic_button = QPushButton("HMI CONT MAP")
        hmic_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        hmic_button.clicked.connect(lambda: self.show_map("conthmi"))
        im_layout.addWidget(hmic_button)
        # Add HMI SAAS Button
        picture_button = QPushButton("SAAS: HMI CONT MAP")
        picture_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        picture_button.clicked.connect(lambda: self.show_map("saashmi"))
        im_layout.addWidget(picture_button)
        # Add buttons to layout
        layout.addLayout(im_layout)
        #----------------------------------------------------------------------
        # Science Table Output
        space_label = QLabel(" ")
        space_label.setAlignment(Qt.AlignCenter)
        space_label.setFixedHeight(20)
        space_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addWidget(space_label)
        table_label = QLabel("SCIENCE AR Target Table")
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
    def transform_target_to_rocket(self,Tx,Ty, obstime,launch_starttime,dtime,**kwargs):
        # Get launch window times        
        foxsi_times = sunpy.time.parse_time(launch_starttime,scale='utc') + \
                      np.linspace(np.floor((lwdur/2)/dtime),lwdur,num=dtime) * u.hour
        # Get observed time of AR
        obs_time   = sunpy.time.parse_time(obstime,scale='utc')
        #----------------------------------------------------------------------
        # Get coordinate of FOXSI at apogee.
        foxsi_loc = astropy.coordinates.EarthLocation(lon=rocket_lon, lat=rocket_lat, height=rocket_altitude + rocket_apogee)
        foxsi_skycoord = SkyCoord(foxsi_loc.get_gcrs(foxsi_times))
        # Find observer of AR
        target_observer = astropy.coordinates.get_body("Earth", time=obs_time)
        #----------------------------------------------------------------------
        # Method of using Stonyhurst Coordinates
        # coord_units = [u.deg,u.deg]
        # frame = HeliographicStonyhurst
        # planning_targets = {"A": SkyCoord(lat=20, lon=-84, unit=coord_units, obstime=obstime, frame=frame),
        #             "B": SkyCoord(lat=26, lon=-5, unit=coord_units, obstime=obstime, frame=frame),
        #             "C": SkyCoord(lat=5, lon=67, unit=coord_units, obstime=obstime, frame=frame),
        #             "D": SkyCoord(lat=-8, lon=-47, unit=coord_units, obstime=obstime, frame=frame),
        #             "E": SkyCoord(lat=8, lon=-2, unit=coord_units, obstime=obstime, frame=frame),
        #             "F": SkyCoord(lat=-11, lon=65, unit=coord_units, obstime=obstime, frame=frame),
        #             "G": SkyCoord(lat=26, lon=-57, unit=coord_units, obstime=obstime, frame=frame),
        #             "H": SkyCoord(lat=11, lon=49, unit=coord_units, obstime=obstime, frame=frame)}
        # planning_targets = list(planning_targets.values())
        #----------------------------------------------------------------------
        # Iterate until AR is on disk in case it is off disk
        new_Tx = np.nan
        factor = 1
        while np.any(np.isnan(new_Tx)):
            # Create Helioprojective coordinates
            planning_target = SkyCoord(Tx * u.arcsec, Ty * u.arcsec,
                              frame=Helioprojective(observer=target_observer))
            #planning_target = planning_targets[cntr]
            # Differentially rotate target to launch time.
            target_diffrot = SkyCoord(\
                             RotatedSunFrame(base=planning_target.frame,\
                                             rotated_time=foxsi_skycoord.obstime,\
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
    def science_target(self,mission_dict):
        #----------------------------------------------------------------------
        self.mission      = mission_dict['mission']
        sparcsoffset = mission_dict['sparcsoffset']
        dtime = 5
        dtimes = np.linspace(np.floor((lwdur/2)/dtime),lwdur,num=dtime)
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
            row_names.append(" "+str(azlist[ri])+". "+parts[0].split('.')[-1]+"         ")
            values = parts[1].replace(' (', ' ').replace(')', '').replace(',', ' ').split()
            fvals  = []
            # Assign HPC Values
            Tx = float(values[-2].replace('"', ''))
            Ty = float(values[-1].replace('"', ''))
            newTx,newTy,launchlon,launchlat = self.transform_target_to_rocket(Tx,Ty,obsdate,launchdate,dtime)
            # Check if on solar disk
            if newTx == Tx:
                fvals.append('Yes')
            else:
                fvals.append('No')
            # Plot in North and West coordinates
            fvals.append("({:+.2f}, {:+.2f})".format(newTy, newTx))
            #------------------------------------------------------------------
            launchcoor = ["({:+.2f}, {:+.2f})".format(x, y) \
                          for x, y, in zip(launchlat.to_value(u.arcsec),launchlon.to_value(u.arcsec))]
            fvals.extend(launchcoor)
            # Assign Column Values in Row
            offsetlaunchcoor = ["({:+.2f}, {:+.2f})".format(x+np.nan_to_num(sparcsoffset[1], nan=0), y+np.nan_to_num(sparcsoffset[0], nan=0)) \
                          for x,y in zip(launchlat.to_value(u.arcsec),launchlon.to_value(u.arcsec))]
            # Add Column Values
            fvals.extend(offsetlaunchcoor)
            table_data.append(fvals)
        #----------------------------------------------------------------------
        # Assign Column Names
        column_names = ["On Disk?",' \nSDO\n(N, W) ["]\nLWT{:+.1f}hr\n   '.format(rot_duration)]
        launchnames  =  [f' \n{self.mission}\nEXP\n(N, W) ["]\nLWT+{i}hr\n   ' for i in dtimes]
        offsetnames  =  [f' \n{self.mission}\nSPARCS\n(N, W) ["]\nLWT+{i}hr\n   ' for i in dtimes]
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
    def sparcs_target(self,mission):
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
        self.new_window = SPARCSGUI(self,mission)
        self.new_window.showMaximized()
    ###########################################################################
    def show_map(self,maptype):
        #----------------------------------------------------------------------
        # Get the input data
        data = self.data_text.toPlainText().strip()
        # Split the data into rows
        rows = data.split('\n')
        #----------------------------------------------------------------------
        # Get AR data
        rows_data = rows[2:]
        # Get AR Observation Date
        obsdate = (rows[0].split(' '))[-1]
        # Get the selected launch date
        launchdate = self.date_time_edit.dateTime().toString("yyyy/MM/ddTHH:mm:ss")
        # Initialize Table contents
        coor_data   = np.zeros((2,len(rows_data)))
        row_names    = []
        #----------------------------------------------------------------------
        # Iterate each AR
        for ri,row in enumerate(rows_data):
            parts = row.split(': ')
            if len(parts) >= 2:
                row_names.append(" "+str(azlist[ri])+". "+parts[0].split('.')[-1]+"         ")
                values = parts[1].replace(' (', ' ').replace(')', '').replace(',', ' ').split()
                # Assign HPC Values
                Tx = float(values[-2].replace('"', ''))
                Ty = float(values[-1].replace('"', ''))
                newTx,newTy,launchlon,launchlat = self.transform_target_to_rocket(Tx,Ty,obsdate,launchdate)
                coor_data[0,ri] = launchlon[2].value
                coor_data[1,ri] = launchlat[2].value

        # Open a file explorer dialog to select a file
        default_dir = os.path.join(os.getcwd(), "Solar_Maps")
        img_dir = os.path.join(default_dir, "Images")
        data_dir = os.path.join(default_dir, "Data")
        # Check if the directory exists, if not create it
        if not os.path.exists(img_dir):
            os.makedirs(img_dir)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        #----------------------------------------------------------------------
        if maptype =='saashmi':
            data_dir = os.path.join(data_dir, "HMI","CONT")
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
            file_path, _ = QFileDialog.getOpenFileName(self, "Select SAAS HMI Data",data_dir, "Image Files (*.png *.jpg *.jpeg *.bmp *.fits)")
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
                # Set the background color of the plot
                ax = fig.add_subplot(projection=hmi_rotated)
                ax.set_facecolor('black')
                ax.imshow(data, cmap='gray',extent=[0,datasz[0],0,datasz[1]])
                # Set Targets on Image
                coor = np.zeros((2,len(rows_data)))
                if row_names is not None:
                    # Get scale
                    scl = hmi_rotated.scale[0].value
                    # Transform to SAAS Frame
                    coor[0,:] = -1*coor_data[1,:]/scl + datasz[1]/2
                    coor[1,:] =  coor_data[0,:]/scl + datasz[0]/2
                
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
                ax.set_title("HMI Continuum"+hmi_map.date.strftime('%Y-%m-%dT%H:%M:%S')+" - SAAS Frame ", fontsize=38)
                # Add annotation at the bottom
                arrow_text = "Solar North"
                ax.annotate(arrow_text, xy=(0.35, -.03), xytext=(0.45, -.04), fontsize=32,
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(facecolor='black', arrowstyle='->', linewidth=4))
                ax.set_xlim([-datasz[0]/3,datasz[0]+datasz[0]/3])
                ax.set_ylim([-datasz[1]/3,datasz[0]+datasz[1]/3])
                if row_names != []:
                    # Plot each circle with a different color
                    for i in range(len(row_names)):
                        plt.scatter(coor[0,i], coor[1,i], marker='o', 
                                    facecolors='none', label=row_names[i],linewidths=5, s=600, color=colors[i])
                    # Add legend
                    plt.legend(fontsize=24)
                plt.show()
                # Save and Open File
                sfile = os.path.join(img_dir, "HMI_IMAGE_"+hmi_map.date.strftime('%Y-%m-%dT%H%M%S'))
                plt.savefig(sfile)
                plt.close()
                subprocess.Popen(sfile+".png", shell=True)
            else:
                # Display the selected picture (placeholder code)
                QMessageBox.warning(self, "Error", "No Image!")
        #----------------------------------------------------------------------
        if maptype =='aia':
            data_dir = os.path.join(data_dir, "AIA")
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
            file_path, _ = QFileDialog.getOpenFileName(self, "Select AIA Data",data_dir, "Image Files (*.png *.jpg *.jpeg *.bmp *.fits)")
            # Check if a file is selected
            if file_path:
                plt.rcParams.update({'font.size':32})
                aiamap = sunpy.map.Map(file_path)
                fig = plt.figure(figsize=(20,20))
                ax1 = fig.add_subplot( projection=aiamap)
                aiamap.plot(axes=ax1, title="AIA MAP "+aiamap.date.strftime('%Y-%m-%dT%H:%M:%S')+" (Solar North)",
                            clip_interval=(60, 99.8)*u.percent)
                if row_names != []:
                    # Plot each circle with a different color
                    for i in range(len(row_names)):
                                
                        ax1.plot_coord(SkyCoord(coor_data[0,i] * u.arcsec, coor_data[1,i] * u.arcsec, frame=aiamap.coordinate_frame), "o",
                                      label=row_names[i],markersize=30,color=colors[i],markerfacecolor='none',markeredgewidth=8)
                                    
                    # Add legend
                    ax1.legend(fontsize=24)
                plt.show()
                # Save and Open File
                sfile = os.path.join(img_dir, "AIA_IMAGE_"+aiamap.date.strftime('%Y-%m-%dT%H%M%S'))
                plt.savefig(sfile)
                plt.close()
                subprocess.Popen(sfile+".png", shell=True)    
            else:
                # Display the selected picture (placeholder code)
                QMessageBox.warning(self, "Error", "No Image!")
        #----------------------------------------------------------------------
        if maptype =='maghmi':
            data_dir = os.path.join(data_dir, "HMI","MAG")
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
            file_path, _ = QFileDialog.getOpenFileName(self, "Select HMI MAG Data",data_dir, "Image Files (*.png *.jpg *.jpeg *.bmp *.fits)")
            # Check if a file is selected
            if file_path:
                plt.rcParams.update({'font.size':32})
                map_hmi = sunpy.map.Map(file_path)
                map_hmi = map_hmi.rotate(order=3)
                map_hmi.plot_settings['cmap'] = "hmimag"
                map_hmi.plot_settings['norm'] = plt.Normalize(-1500, 1500)
                fig = plt.figure(figsize=(20,20))
                ax1 = fig.add_subplot( projection=map_hmi)
                map_hmi.plot(axes=ax1, title="HMI MAGNETOGRAM MAP "+map_hmi.date.strftime('%Y-%m-%dT%H:%M:%S')+" (Solar North)")
                            #clip_interval=(10, 99.99)*u.percent)
                if row_names != []:
                    # Plot each circle with a different color
                    for i in range(len(row_names)):
                                
                        ax1.plot_coord(SkyCoord(coor_data[0,i] * u.arcsec, coor_data[1,i] * u.arcsec, frame=map_hmi.coordinate_frame), "o",
                                      label=row_names[i],markersize=30,color=colors[i],markerfacecolor='none',markeredgewidth=8)
                                    
                    # Add legend
                    ax1.legend(fontsize=24)
                plt.show()
                # Save and Open File
                sfile = os.path.join(img_dir, "HMI_MAGNETIC_IMAGE_"+map_hmi.date.strftime('%Y-%m-%dT%H%M%S'))
                plt.savefig(sfile)
                plt.close()
                subprocess.Popen(sfile+".png", shell=True) 
            else:
                # Display the selected picture (placeholder code)
                QMessageBox.warning(self, "Error", "No Image!")
        #----------------------------------------------------------------------
        if maptype =='conthmi':
            data_dir = os.path.join(data_dir, "HMI","CONT")
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
            file_path, _ = QFileDialog.getOpenFileName(self, "Select HMI CONT Data",data_dir, "Image Files (*.png *.jpg *.jpeg *.bmp *.fits)")
            # Check if a file is selected
            if file_path:
                plt.rcParams.update({'font.size':32})
                map_hmi = sunpy.map.Map(file_path)
                map_hmi = map_hmi.rotate(order=3)

                #map_hmi.plot_settings['cmap'] = "hmimag"
                #map_hmi.plot_settings['norm'] = plt.Normalize(-1500, 1500)
                fig = plt.figure(figsize=(20,20))
                ax1 = fig.add_subplot( projection=map_hmi)
                map_hmi.plot(axes=ax1, title="HMI Continuum MAP "+map_hmi.date.strftime('%Y-%m-%dT%H:%M:%S')+" (Solar North)")
                            #clip_interval=(10, 99.99)*u.percent)
                
                if row_names != []:
                    # Plot each circle with a different color
                    for i in range(len(row_names)):
                                
                        ax1.plot_coord(SkyCoord(coor_data[0,i] * u.arcsec, coor_data[1,i] * u.arcsec, frame=map_hmi.coordinate_frame), "o",
                                      label=row_names[i],markersize=30,color=colors[i],markerfacecolor='none',markeredgewidth=8)
                                    
                    # Add legend
                    ax1.legend(fontsize=24)
                plt.show()
                # Save and Open File
                sfile = os.path.join(img_dir, "HMI_CONT_IMAGE_"+map_hmi.date.strftime('%Y-%m-%dT%H%M%S'))
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
        filename = os.path.join(filedir, self.mission+"-Science-Target_Tables_"+launchdate+".csv")
        # Access the table model
        model = self.table_view.model()
        if model is not None:
            # Open the CSV file for writing
            with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                # Write column headers
                headers = ['Target']
                headers.extend([model.headerData(col, Qt.Horizontal).replace("OPTICS","FOXSI "
                                                                    ).replace("EXP"," EXP"
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
    def __init__(self, parent,mission_dict):
        super().__init__()
        self.setWindowTitle(mission_dict['mission'] +" SPARCS Targeting")
        self.initUI(parent,mission_dict)
    ###########################################################################
    def initUI(self,parent,mission_dict):
        self.mission      = mission_dict['mission']
        sparcsoffset = mission_dict['sparcsoffset']
        dtime        = mission_dict['dtime']
        dtimes = np.linspace(np.floor((lwdur/2)/dtime),lwdur,num=dtime)
        midt = int(len(dtimes)/2)
        dlon = mission_dict['dlon']
        dlat = mission_dict['dlat']
        # Create Layout
        layout = QVBoxLayout()
        # Horizontal layout for text edit and date input
        hbox = QHBoxLayout()
        title_label = QLabel(self.mission +" SPARCS Targeting")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addLayout(hbox)
        # Add title label
        title_label = QLabel(self.mission +" SPARCS Targeting")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 48px; font-weight: bold;")
        layout.addWidget(title_label)
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
            row_names.append(" "+str(azlist[ri])+". "+parts[0].split('.')[-1]+"         ")
            values = parts[1].replace(' (', ' ').replace(')', '').replace(',', ' ').split()
            fvals  = []
            # Assign HPC Values
            Tx = float(values[-2].replace('"', ''))
            Ty = float(values[-1].replace('"', ''))
            newTx,newTy,launchlon,launchlat = parent.transform_target_to_rocket(Tx,Ty,obsdate,launchdate,dtime)
            #------------------------------------------------------------------
            # Set center values
            fvals.extend(["(  {:+.2f}, {:+.2f}  )".format(launchlat.to_value(u.arcsec)[midt],launchlon.to_value(u.arcsec)[midt])])
            fvals.extend(["(  {:+.2f}, {:+.2f}  )".format(sparcsoffset[1],sparcsoffset[0])])
            # Compute new delta values
            lonround = np.round(launchlon.to_value(u.arcsec)+np.nan_to_num(sparcsoffset[0], nan=0))
            latround = np.round(launchlat.to_value(u.arcsec)+np.nan_to_num(sparcsoffset[1], nan=0))
            lonsprc = [np.round((x - lonround[midt]) / dlon) * dlon for x in lonround]
            latsprc = [np.round((x - latround[midt]) / dlat) * dlat for x in latround]
            lonsprc[midt] = lonround[midt]
            latsprc[midt] = latround[midt]
            # Assign Column Values in Row
            offsetlaunchcoor = ["{:+}".format(x) for x in lonsprc]
            offsetlaunchcoor[midt] = "(  {:+.2f}, {:+.2f}  )".format(latround[midt],lonround[midt])
            # Add Column Values
            fvals.extend(offsetlaunchcoor)
            fvals.extend(['East'])
            table_data.append(fvals)
        #----------------------------------------------------------------------
        # Assign Column Names
        column_names = [f' \n{self.mission}\nEXP\n(N, W) ["]\nLWT+2hr\n   ']
        column_names.extend([f' \n{self.mission}\nEXP-SPARCS Offset\n(ΔN, ΔW) ["]\n   '])
        
        column_names.extend([f' \n{self.mission}\nSPARCS\nΔW ["]\nLWT+{i}hr\n   ' for i in dtimes[0:(midt)]])
        column_names.extend([f' \n{self.mission}\nSPARCS\n( N, W ) ["]\nLWT+2hr\n   '])
        column_names.extend([f' \n{self.mission}\nSPARCS\nΔW ["]\nLWT+{i}hr\n   ' for i in dtimes[(midt+1):]])
        column_names.extend(['Payload\nOrientation'])
        # Save payload orientation
        self.payload_ornt = ['east']*len(row_names)
        # Check for Errors
        if not row_names or not column_names or not table_data:
            QMessageBox.warning(self, "Error", "Data could not be processed.")
            return
        #----------------------------------------------------------------------
        # Set up table
        self.model = QStandardItemModel(len(row_names), len(column_names))
        self.model.setHorizontalHeaderLabels(column_names)
        self.model.setVerticalHeaderLabels(row_names)
        for row_idx, row in enumerate(table_data):
            # Extract every other value as it corresponds to the data we need
            for col_idx, value in enumerate(row):
                item = QStandardItem(value)
                item.setTextAlignment(Qt.AlignCenter)
                self.model.setItem(row_idx, col_idx, item)
        # Set font style for row names and column names
        font = QFont()
        font.setBold(True)
        for row_idx in range(len(row_names)):
            self.model.setHeaderData(row_idx, Qt.Vertical, font, role=Qt.FontRole)
        for col_idx in range(len(column_names)):
            self.model.setHeaderData(col_idx, Qt.Horizontal, font, role=Qt.FontRole)
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
        self.table_view.setModel(self.model)   
        self.model.dataChanged.connect(self.handleClicked)
        #self.setCentralWidget(self.table_view)
        #----------------------------------------------------------------------
        # Export to CSV
        convert_button = QPushButton("Export to CSV")
        convert_button.setStyleSheet("font-size: 38px; font-weight: bold;")
        convert_button.clicked.connect(lambda: self.sparcs_exportToCSV(parent.date_time_edit.dateTime(),dtimes))
        layout.addWidget(convert_button)
    ###########################################################################
    def handleClicked(self, index):
        model = self.table_view.model()
        colnum = model.columnCount()
        if index.column() == (colnum-1):  # Check if clicked item is in the first column
            edited_row = index.row()
            # Get new value
            edited_value = (model.item(edited_row,model.columnCount()-1).text()) 
            if edited_value.lower() == 'east':
                # Change values
                delcol = colnum-4
                for dc in range(delcol):
                    # Update corresponding cell in the second column
                    mult = -1
                    if dc >= delcol/2:
                        dc = dc+1
                        mult = 1
                    new_value = mult*abs(float((model.item(edited_row,dc+2).text())))
                    new_item = QStandardItem("{:+}".format(new_value))
                    new_item.setTextAlignment(Qt.AlignCenter)
                    model.setItem(edited_row, dc+2, new_item)
            if edited_value.lower() == 'west':
                # Change values
                delcol = colnum-4
                for dc in range(delcol):
                    # Update corresponding cell in the second column
                    mult = 1
                    if dc >= delcol/2:
                        dc = dc+1
                        mult = -1
                    new_value = mult*abs(float((model.item(edited_row,dc+2).text())))
                    new_item = QStandardItem("{:+}".format(new_value))
                    new_item.setTextAlignment(Qt.AlignCenter)
                    model.setItem(edited_row, dc+2, new_item)

    ###########################################################################
    def sparcs_exportToCSV(self,launchdate,dtimes):
        launchtime = launchdate.toString("yyyy/MM/ddTHH:mm:ss")
        launchday  = launchdate.toString("yyyy_MM_dd")
        # Create filename
        filedir  = os.path.join(os.getcwd(), "SPARCS-Target_Tables")
        # Check if the directory exists, if not create it
        if not os.path.exists(filedir):
            os.makedirs(filedir)
        filename = os.path.join(filedir, self.mission+"-SPARCS-Target_Tables_"+launchday+".csv")
        # Access the table model
        model = self.table_view.model()
        # Open the CSV file for writing
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Add Header
            headers = [launchdate.toString("yyyy-MM-dd")]
            headers.extend([model.headerData(col, Qt.Horizontal).replace("\n"," "
                                                               ).replace("Δ"," d"
                                                               ).replace("hr","hr)"
                                                               ).replace("LWT"," (LWT") for col in range(model.columnCount())])
            writer.writerow(headers)
            
            
            # Add UTC Time
            time_utcdata = ["Start Time [UTC]"," "," "]
            time_utc     = sunpy.time.parse_time(launchtime,scale='utc') + dtimes * u.hour
            time_utc_arr = [f"{t.hour:02d}:{t.minute:02d}" for t in time_utc.to_datetime()]
            time_utcdata.extend(time_utc_arr)
            time_utcdata.extend([" "])
            writer.writerow(time_utcdata)
            # Add Local Time
            time_lcldata = ["Start Time [AKDT]"," "," "]
            time_lcl     = sunpy.time.parse_time(launchtime,scale='utc') + (dtimes+lcl_time) * u.hour
            time_lcl_arr = [f"{t.hour:02d}:{t.minute:02d}" for t in time_lcl.to_datetime()]
            time_lcldata.extend(time_lcl_arr)
            time_lcldata.extend([""])
            writer.writerow(time_lcldata)
            #Spacer
            writer.writerow([" "])
            # Write column headers
            headers = ['Target']
            writer.writerow(headers)
            # Write data rows
            for row in range(model.rowCount()):
                row_name = model.headerData(row, Qt.Vertical)
                row_data = [row_name] + [model.index(row, col).data() for col in range(model.columnCount())]
                writer.writerow(row_data)
        QMessageBox.information(self, 'CSV Generated', 'CSV File Created!')
###############################################################################
def create_application():
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app
###############################################################################
###############################################################################
###############################################################################
# Main 
if __name__ == '__main__':
    app = create_application()
    converter = FOXSITargetGUI()
    # Maximize the window
    converter.showMaximized()  
    sys.exit(app.exec_())