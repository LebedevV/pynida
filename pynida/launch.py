#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import scipy
from scipy import ndimage#, interpolate
from scipy.optimize import fsolve

import pandas as pd
#for img rotation
#from PIL import Image

from matplotlib.widgets import RectangleSelector
#import matplotlib.image as mpimg
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

#For areas segmentation only
from skimage import measure #<0.16
#

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QApplication, QFileDialog, QMessageBox #, QLineEdit, QWidget, QInputDialog,  QMainWindow
from PyQt5.QtGui import QPixmap, QDoubleValidator, QIntValidator #QIcon,

import cv2

# Import of the graphical template
from pynida.frontend import Ui_MainWindow
from pynida.backend import oneshot_contours, oneshot_areas_select, oneshot_areas_proc
from pynida.simple_functions import linear, sneddon, Convert, eab, parse_corr_data, parse_hys_data, get_n_frame, parse_video_data, parse_old_data
from pynida.frames_preview import frames_preview

from pynida.window_comparison import *

# Main window
class mywindow(QtWidgets.QMainWindow):
	def __init__(self):
		super(mywindow, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self)

		units = ['nm','μm','mm','m']
		for i in units:
			self.ui.comboBox.addItem(i)
		#self.ui.radioButton.setChecked(True)
		global progress
		global progress_area
		progress = self.ui.progressBar
		progress_area = self.ui.progressBar_2

		self.separator = '---' # for video data parsing
		self.parameters = None
		self.savedir = None
		self.name_of_videofile = ''

		#limits user input for float values. Top-bottom values are not checking
		for i in ['startval1','startval2','stepcoeff1','stepcoeff2',
				'AreaScale','AreaScale_2','fps','shift','pnm2','8']:
			getattr(self.ui,'lineEdit_' + i).setValidator(QDoubleValidator())
		for i in ['force_scale','k_tip','tip_radius_2','tip_radius',
				'sphere_radius','cyl_radius','probe_modulus','substr_modulus',
				'probe_poiss','substr_poiss','RotationAngle','lineEdit']:
			getattr(self.ui,i).setValidator(QDoubleValidator())

		#limits user input for int values. Top-bottom values are not checking
		for i in ['11','10','20','9','2','4','5','TestFrN','minarea_1','minarea_1',
				'14', '15', '17','18']:#'TestFrN_2'
			getattr(self.ui,'lineEdit_' + i).setValidator(QIntValidator())


		#click-slots connections for all buttons
		self.ui.pushButton_start_tracker.clicked.connect(self.btnClicked_start_tracker)
		self.ui.pushButton_3.clicked.connect(self.btnClicked_P3)
		self.ui.pushButton_7.clicked.connect(self.btnClicked_P7)
		self.ui.pushButton_4.clicked.connect(self.btnClicked_P4)
		self.ui.pushButton_SelectScalebar.clicked.connect(self.btnClicked_SelectScalebar)
		self.ui.pushButton_5.clicked.connect(self.btnClicked_show)

		self.ui.pushButton_6.clicked.connect(self.btnClicked_selectfile)
		self.ui.pushButton_11.clicked.connect(self.btnClicked_selectfile)
		self.ui.pushButton_12.clicked.connect(self.btnClicked_selectfile)

		self.ui.pushButton_9.clicked.connect(self.btnClicked_select_saveDir)
		self.ui.pushButton_14.clicked.connect(self.btnClicked_select_saveDir)
		self.ui.pushButton_10.clicked.connect(self.btnClicked_select_saveDir)
		self.ui.pushButton_13.clicked.connect(self.btnClicked_select_saveDir)

		self.ui.pushButton_8.clicked.connect(self.btnClicked_frames)

		self.ui.checkBox_3.clicked.connect(self.checkBox_allow_out)
		self.ui.checkBox_7.clicked.connect(self.checkBox_allow_out)
		self.ui.checkBox_8.clicked.connect(self.checkBox_allow_out)

		#self.ui.checkBox_Rotate.clicked.connect(self.allow_dial)

		self.ui.pushButton_crop.clicked.connect(self.crop_dialog)
		self.ui.pushButton_crop_2.clicked.connect(self.crop_dialog)
		self.ui.pushButton_crop_3.clicked.connect(self.crop_dialog)

		self.ui.button_hys_data_import.clicked.connect(self.btnClicked_import_hysdata)
		self.ui.button_video_data_import.clicked.connect(self.btnClicked_import_videodata)
		self.ui.button_corr_data_import.clicked.connect(self.btnClicked_import_corrdata)
		self.ui.button_old_data_import.clicked.connect(self.btnClicked_import_olddata)

		self.ui.radioButton_pnmbybar.clicked.connect(self.radiobutton_check)
		self.ui.radioButton_pnmbymag.clicked.connect(self.radiobutton_check)

		self.ui.radioButton_Y_lin.clicked.connect(self.radiobutton_check_Young)
		self.ui.radioButton_Y_Cyl.clicked.connect(self.radiobutton_check_Young)
		self.ui.radioButton_Y_Sph.clicked.connect(self.radiobutton_check_Young)

#		self.ui.checkBox_Cont.clicked.connect(self.checkBox_ArC)
#		self.ui.checkBox_Area.clicked.connect(self.checkBox_ArC)
		self.ui.pushButton_TestCall.clicked.connect(self.btnClicked_TestCall)
		self.ui.pushButton_TestCall_2.clicked.connect(self.btnClicked_TestCall_contours)
		self.ui.pushButton_StartSegm.clicked.connect(self.btnClicked_FullAreaAnalysis)
		self.ui.pushButton_StartSegm_2.clicked.connect(self.btnClicked_FullContoursAnalysis)
		self.ui.pushButton_AreasOpen.clicked.connect(self.btnClicked_AreasOpen)
		self.ui.pushButton_ContoursOpen.clicked.connect(self.btnClicked_ContoursOpen)
		self.ui.pushButton_ConvertAreas.clicked.connect(self.btnClicked_AreasConv)
		self.ui.pushButton_ConvertContours.clicked.connect(self.btnClicked_ContoursConv)
		self.ui.pushButton_FinPlot.clicked.connect(self.plot_converted)
		self.ui.pushButton_15.clicked.connect(self.plot_old)

		self.ui.pushButton_openComparison.clicked.connect(self.btnClicked_openComparison)
		self.ui.checkBox_Rotate.toggled.connect(self.areas_recall_part)
		self.ui.RotationAngle.textChanged.connect(self.areas_recall_part)
		self.ui.lineEdit_9.textChanged.connect(self.areas_recall_part)
		

	def areas_recall_part(self):
		self.ui.label_10.clear()
		self.ui.label_13.clear()
		self.ui.label_10.setText('Tip movement')
		self.ui.label_13.setText('Sample movement')
		progress.setValue(0)
		for i in ['pushButton_start_tracker']:
			getattr(self.ui,i).setEnabled(False)

#		if self.ui.checkBox_Rotate.isChecked():
		todo = True
		try:
			angle = float(self.ui.RotationAngle.text())
		except:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
		if todo:
			try:
				cap = cv2.VideoCapture(self.name_of_videofile)
				len_cap = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
				i = 0
				while i < len_cap:
					try:
						_,img = cap.read()
						break
					except:
						print(i)
					i+=1
				cap.release()

				if self.ui.checkBox_Rotate.isChecked():
					img = ndimage.rotate(img, angle)

				self.update_crop(img)
			except:
				print("Some error with crop area")

	def areas_recall(self):
		self.ui.label_10.clear()
		self.ui.label_13.clear()
		self.ui.label_10.setText('Tip movement')
		self.ui.label_13.setText('Sample movement')
		self.ui.radioButton_pnmbymag.setChecked(True)
		progress.setValue(0)
		for i in ['pushButton_start_tracker','radioButton_pnmbybar','lineEdit_2','comboBox']:
			getattr(self.ui,i).setEnabled(False)


	def btnClicked_frames(self):

		a = self.ui.lineEdit_11.text()
		b = self.ui.lineEdit_10.text()
		ln = self.ui.lineEdit_6.text()
		angle = self.ui.RotationAngle.text()
		todo = False
		try:
			a, b, ln = int(a), int(b), int(ln)
			if a > b:
				a,b = b,a

			if a < 0: a = 0
			if b < 0: b = 0
			if a > ln: a = ln
			if b > ln: b = ln

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame numbers")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
			todo = False

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo and abs(a-b) > 100:
			qm = QMessageBox
			ret = qm.question(self,'', "You are about to load more than 100 frames simultaneously. \n Are you sure?", qm.Yes | qm.No)
			if ret == qm.No:
				todo = False

		if todo:
			try:
				frames_preview(self.name_of_videofile,
						a, b, ln,
						self.ui.checkBox_Rotate.isChecked(), angle,
						self.ui.checkBox_view_crop.isChecked(), (crop_x0, crop_x1, crop_y0, crop_y1))
			except:
				QMessageBox.about(self, "Unexpected error", "Please, report this case")


	def plot_converted(self):
		plt.figure(figsize = (8, 6))
		grid = plt.GridSpec(22, 4, hspace = 0., wspace = 0.)
		ax0=plt.subplot(grid[0:5, 0:])

		ax00 = plt.subplot(grid[5:10, 0:], sharex = ax0)
		ax01 = ax00.twinx()
		ax00.get_yaxis().set_visible(False)
		ax02 = plt.subplot(grid[10:15, 0:], sharex = ax0)
		for tt in [ax0,ax02,ax01,ax00]:
			for x in tt.get_xticklabels():
				x.set_visible(False)
			for x in tt.get_xticklines():
				x.set_visible(False)
		ax1 = plt.subplot(grid[16:22, 0:], sharex = ax0)
		ax2 = ax1.twinx()
		#ax4=ax0.twinx()
		ax02.set_ylabel("$Volume, \\mu m^3$", color = 'r')
		#ax0.set_xlabel("$t$, $s$")

		ax0.set_ylabel("$Width$, $nm$", color = 'b')
		ax01.set_ylabel("$Height$, $nm$", color= 'g')
		#ax02.set_ylabel("$Len$, $nm$", color= 'b')
		ax1.set_xlabel("$t$, $s$")
		ax2.set_ylabel("$h$, $nm$", color= 'r')
		ax1.set_ylabel("$F$, $\mu N$", color= 'b')
		plt.subplots_adjust(right=0.85, left=0.15, top=0.9, bottom=0.13,hspace=0,wspace=0)

		for tt in [ax0,ax00,ax01,ax02,ax1,ax2]:
			tt.grid(True)
			tt.xaxis.label.set_size(16)
			tt.yaxis.label.set_size(16)
			tt.tick_params(labelsize=14)
		ax2.grid(False)


#		ax0.plot(ntsh,nmid, "ob", label= 'Width (mid)',marker= '.')
		ax0.plot(self.a_frame_time,self.a_sqy, "ob", label= 'Width (box)',marker= '.')

	#	ax0.plot(time_sh,boxx, "-og", label= 'F(d)')
		#ax01.plot(time_sh,boxx, "og", label= 'Height',marker= '.')
		if self.ui.lineEdit_29.text():
			ax01.plot(self.a_frame_time,self.a_sqx, "og", label= 'Height (box)',marker= '.')
			ax02.plot(self.c_frame_time,self.c_vol, "or", label= 'Volume',marker= '.')
		#ax02.plot(ntsh,nps, "ob", label= 'Volume',marker= '.')

		ax1.plot(self.fin_ht,self.fin_hf,"b")
		ax2.plot(self.fin_ht,self.fin_hx,"r")
		plt.savefig(self.ui.lineEdit_27.text()[:-5]+'_visall.png')

#to update: move to GUI
		ax0.set_xlim(150,200)
#
		plt.savefig(self.ui.lineEdit_27.text()[:-5]+'_visall_cut_.png')
		plt.close()
		self.ui.label_finimg.setPixmap(QPixmap(self.ui.lineEdit_27.text()[:-5] + '_visall.png').scaled(self.ui.label_finimg.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
		self.ui.label_finimg.show()

	def btnClicked_AreasConv(self):
		df = pd.read_csv(self.area_data_import)
		try:
			header_file_name = self.area_data_import.split('.')[0] + '.shpx'
			with open(header_file_name, 'r') as file:
				parameters = eval(file.read())
			a_scale = parameters['area_scale']
		except:
			a_scale = 1

		areas_conv = False
		try:
			tmp_multiplier = float(self.ui.lineEdit_pnm2.text()) / a_scale
			self.a_frame_time = (df['frame_number'].values + float(self.ui.lineEdit_shift.text())) * float(self.ui.lineEdit_fps.text())
			self.a_sqx = df['square_x'].values * tmp_multiplier
			self.a_sqy = df['square_y'].values * tmp_multiplier
			df['tip_sample_contact_length'] = df['tip_sample_contact_length'] * tmp_multiplier
			df['sample_surface_contact_length'] = df['sample_surface_contact_length'] * tmp_multiplier
			df['mid_projection_x'] = df['mid_projection_x'] * tmp_multiplier
			df['mid_projection_y'] = df['mid_projection_y'] * tmp_multiplier
			df['area'] = df['area'] * tmp_multiplier ** 2
			areas_conv = True
		except:
			print("Error in the area data conversion")

		if self.ui.pushButton_ConvertContours.isChecked() and areas_conv:
			self.ui.pushButton_FinPlot.setEnabled(True)

		if areas_conv:
			plt.plot(self.a_frame_time, self.a_sqy, label='X size of area')
			plt.plot(self.a_frame_time, self.a_sqx, label='Y size of area')
			plt.ylabel('d, nm')
			plt.xlabel('t, s')
			plt.legend()
			plt.savefig(self.ui.lineEdit_26.text()[:-5] + '_XYbox_sh.png')
			plt.close()

			plt.plot(self.a_frame_time, df['tip_sample_contact_length'], label='Length of Tip - Sample contact')
			plt.plot(self.a_frame_time, df['sample_surface_contact_length'], label='Length of Sample - Surface contact')
			plt.ylabel('d, nm')
			plt.xlabel('t, s')
			plt.legend()
			plt.savefig(self.ui.lineEdit_26.text()[:-5] + '_lts_sh.png')
			plt.close()

			plt.plot(self.a_frame_time, df['mid_projection_x'], label='Midline X')
			plt.plot(self.a_frame_time, df['mid_projection_y'], label='Midline Y')
			plt.ylabel('d, nm')
			plt.xlabel('t, s')
			plt.legend()
			plt.savefig(self.ui.lineEdit_26.text()[:-5] + '_mid_sh.png')
			plt.close()

			plt.plot(self.a_frame_time, df['area'], label='Midline X')
			plt.ylabel('S, nm$^2$')
			plt.xlabel('t, s')
			plt.legend()
			plt.savefig(self.ui.lineEdit_26.text()[:-5] + '_A_sh.png')
			plt.close()

	def btnClicked_ContoursConv(self):
		d = open(self.contours_data_import,"r")
		c_scale = 1
		while 1:
			s = d.readline()
			if s.startswith(self.separator):
				break
			if 'scale' in s:
				s=s.split('\t')
				c_scale = float(s[1])
		c_data_desc = d.readline()
		c_data = d.readlines()
		d.close()
		c_data_v = [s.split('\t') for s in c_data]
		c_data = np.array(c_data_v,dtype= 'float')
		#print(a_data[:,0])

		try:
#		if 1:
			c_data_desc=c_data_desc.split('\t')
			i = 0
			while i<len(c_data_desc):
				if 'Frame' in c_data_desc[i]:
					c_frame = c_data[:,i]
				if 'Vol' in c_data_desc[i] and 'top' in c_data_desc[i]:
					c_vol_top = c_data[:,i]
				#if 'corrvol' in c_data_desc[i] and 'top' in c_data_desc[i]:
				#	c_corrvol_top = c_data[:,i]
				#if 'corrvol' in c_data_desc[i] and 'b' in c_data_desc[i]:
				#	c_corrvol_bot = c_data[:,i]
				if 'Vol' in c_data_desc[i] and 'bot' in c_data_desc[i]:
					c_vol_bot = c_data[:,i]
					#self.a_sqy = a_sqy/a_scale*float(self.ui.lineEdit_pnm2.text())
				i+=1
			self.c_vol = (c_vol_top+c_vol_bot)/2./c_scale**3*float(self.ui.lineEdit_pnm2.text())**3/10**9 #vol in mkm^3
			#removal of fake volumes
			self.c_vol[self.c_vol == 0.] = np.nan
			self.c_frame_time = (c_frame +float(self.ui.lineEdit_shift.text()))*float(self.ui.lineEdit_fps.text())

			'''
			plt.close()
			plt.plot(self.c_frame_time,c_vol_top,label= 'Top volume',color= 'green')
			plt.plot(self.c_frame_time,c_vol_bot,label= 'Bot volume',color= 'blue')
			plt.plot(self.c_frame_time,c_vol_bot/2.+c_vol_top/2.,label= 'Avg volume',color= 'red')
			#plt.plot(fr,np.array(D_sqa)*np.pi*np.array(D_sqb)**2,color= 'red')
			plt.legend()
			plt.savefig(self.ui.lineEdit_27.text()[:-5]+'_volume.png')
			plt.close()

			plt.close()
			plt.plot(self.c_frame_time,c_corrvol_top,label= 'Top volume',color= 'green')
			plt.plot(self.c_frame_time,c_corrvol_bot,label= 'Bot volume',color= 'blue')
			plt.plot(self.c_frame_time,c_corrvol_bot/2.+c_corrvol_top/2.,label= 'Avg volume',color= 'red')
			#plt.plot(fr,np.array(D_sqa)*np.pi*np.array(D_sqb)**2,color= 'red')
			plt.legend()
			plt.savefig(self.ui.lineEdit_27.text()[:-5]+'_cvolume.png')
			plt.close()

			plt.close()
			plt.plot(self.c_frame_time,c_vol_top+c_corrvol_top,label= 'Top volume',color= 'green')
			plt.plot(self.c_frame_time,c_vol_bot+c_corrvol_bot,label= 'Bot volume',color= 'blue')
			plt.plot(self.c_frame_time,c_vol_bot/2.+c_vol_top/2.+c_corrvol_bot/2.+c_corrvol_top/2.,label= 'Avg volume',color= 'red')
			#plt.plot(fr,np.array(D_sqa)*np.pi*np.array(D_sqb)**2,color= 'red')
			plt.legend()
			plt.savefig(self.ui.lineEdit_27.text()[:-5]+'_corr_volume.png')
			plt.close()
#			'''
		except:
			print("Error in the contours data conversion")

		if self.ui.pushButton_ConvertAreas.isChecked():
			#try:
			self.ui.pushButton_FinPlot.setEnabled(True)




	def btnClicked_SelectScalebar(self):
		plt.close()

		a = self.ui.lineEdit_9.text()
		ln = self.ui.lineEdit_6.text()
		todo = False
		try:
			a, ln = int(a), int(ln)
			if a < 0: a = 0
			if a > ln: a = ln
			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame number")
			todo = False

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo:
			try:
				global img
				img = get_n_frame(self.name_of_videofile, 0, ln, a)
				img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

				self.ui.radioButton_pnmbybar.setEnabled(True)
				global fig
				global current_ax_scalebar
				fig, current_ax_scalebar = plt.subplots()
				plt.imshow(img.copy(), cmap= 'gray')
				#global rectangle_coordinates
				#rectangle_coordinates = []
				toggle_selector_bar.RS = RectangleSelector(current_ax_scalebar, line_select_callback_scbar,
								   drawtype= 'box', useblit=True,
								   button= [1, 3],  # don't use middle button
								   minspanx=5, minspany=5,
								   spancoords= 'pixels',
								   interactive=True)
				current_ax_scalebar.set_title('Select the scalebar and press "Enter"')
				plt.connect('key_press_event', toggle_selector_bar)
				current_ax_scalebar.axis('off')
				plt.show()
				#print("test")
				#QMessageBox.about(self, "Alarm!", "Select a video file!")

			except:
				QMessageBox.about(self, "Image error", "Please, try another frame")

	def crop_dialog(self):
		plt.close()

		a = self.ui.lineEdit_20.text()
		ln = self.ui.lineEdit_6.text()
		angle = self.ui.RotationAngle.text()
		todo = False
		try:
			a, ln = int(a), int(ln)
			if a < 0: a = 0
			if a > ln: a = ln
			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame number")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if self.ui.lineEdit_3.text():
			try:
				img3 = get_n_frame(self.name_of_videofile,0,ln,a)
				img3 = cv2.cvtColor(img3, cv2.COLOR_BGR2GRAY)

				self.ui.checkBox_view_crop.setEnabled(True)
				self.ui.checkBox_select_crop.setEnabled(True)
				self.ui.checkBox_select_crop.setChecked(True)

				if self.ui.checkBox_Rotate.isChecked():
					img3 = ndimage.rotate(img3, float(self.ui.RotationAngle.text()))
					self.update_crop(img3)
				global fig3
				global current_ax3
				fig3, current_ax3 = plt.subplots()
				plt.imshow(img3.copy(), cmap= 'gray')
				#global rectangle_coordinates
				#rectangle_coordinates = []
				toggle_selector_crop.RS = RectangleSelector(current_ax3, line_select_callback_crop,
												   drawtype= 'box', useblit=True,
												   button= [1, 3],  # don't use middle button
												   minspanx=5, minspany=5,
												   spancoords= 'pixels',
												   interactive=True)
				current_ax3.set_title('Select the area for crop and press "Enter"')
				plt.connect('key_press_event', toggle_selector_crop)
				current_ax3.axis('off')
				plt.show()
			except:
				QMessageBox.about(self, "Error", "Please, try another frame")
		else:
			QMessageBox.about(self, "Alarm!", "Select a video file!")


#	def allow_dial(self):
#		if self.ui.checkBox_Rotate.isChecked():
#			self.ui.dial.setEnabled(True)
#		else:
#			self.ui.dial.setValue(0)
#			self.ui.dial.setEnabled(False)

	def btnClicked_FullAreaAnalysis(self):
		plt.close()

		a = self.ui.lineEdit_14.text()
		b = self.ui.lineEdit_15.text()
		ln = self.ui.lineEdit_12.text()
		angle = self.ui.RotationAngle.text()

		todo = False
		try:
			a, b, ln = int(a), int(b), int(ln)
			if a < 0: a = 0
			if a > ln: a = ln

			if a > b: a, b = b, a
			if a < 0: a = 0
			if b < 0: b = 0
			if b > ln: b = ln
			if a > ln: a = ln
				#it was ln - 1 initially!

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame number")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
			todo = False

		try:
			start_coeff2=float(self.ui.lineEdit_startval2.text())
			min_area2=int(self.ui.lineEdit_minarea_2.text())
			step_coeff2=float(self.ui.lineEdit_stepcoeff2.text())
			sm_blur_val = int(self.ui.lineEdit_sm_blur_val.text())
			start_coeff1=float(self.ui.lineEdit_startval1.text())
			min_area1=int(self.ui.lineEdit_minarea_1.text())
			step_coeff1=float(self.ui.lineEdit_stepcoeff1.text())
			scale = float(self.ui.lineEdit_AreaScale.text())
		except:
			QMessageBox.about(self, "Input error", "Please, check parameters")
			todo = False

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo:
			try:
				for q in [self.ui.tab1,self.ui.tab3,self.ui.tab4,self.ui.tab_3,self.ui.lineEdit_AreaScale,
						self.ui.lineEdit_14,self.ui.lineEdit_15,
						self.ui.lineEdit_TestFrN,self.ui.lineEdit_sm_blur_val,self.ui.pushButton_10,self.ui.pushButton_11,
						self.ui.lineEdit_minarea_1,self.ui.lineEdit_startval1,self.ui.lineEdit_stepcoeff1,
						self.ui.lineEdit_minarea_2,self.ui.lineEdit_startval2,self.ui.lineEdit_stepcoeff2,
						#self.ui.horizontalSlider,
						self.ui.pushButton_TestCall,self.ui.pushButton_StartSegm,
						self.ui.pushButton_a1_b,self.ui.pushButton_a1_t,self.ui.pushButton_a1_l,self.ui.pushButton_a1_r,
						self.ui.pushButton_a2_b,self.ui.pushButton_a2_t,self.ui.pushButton_a2_l,self.ui.pushButton_a2_r,
						self.ui.pushButton_crop_2, self.ui.lineEdit_30, self.ui.checkBox_Rotate_2, self.ui.RotationAngle_2]:
					q.setEnabled(False)

			#create global boolean variables for two sets of "Preserve borders" buttons
				borders_remove1 = []
				borders_remove2 = []
				for i in ['l','r','b','t']:
					borders_remove1.append(getattr(self.ui,"pushButton_a1_"+i).isChecked())
					borders_remove2.append(getattr(self.ui,"pushButton_a2_"+i).isChecked())

				first_fr_num = a
				last_fr_num = b
				cap = cv2.VideoCapture(self.name_of_videofile)
				# len_cap = ln

				self.ui.pushButton_BreakSegm.setEnabled(True)
				self.ui.pushButton_BreakSegm.setChecked(False)

				with open(f"{self.savepath}_areas.shpx",'w') as area_output:
					parameters = {}
					if self.ui.checkBox_Rotate.isChecked():
						parameters['rotation'] = float(self.ui.RotationAngle.text())
					parameters['area_scale'] = float(self.ui.lineEdit_AreaScale.text())
					area_output.write(str(parameters))

				df_data_areas = pd.DataFrame()
				i = 0
				while i <= last_fr_num:
					self.ui.lineEdit_CurrentFrameAreas.setText(str(i))
					QApplication.processEvents()
					if self.ui.pushButton_BreakSegm.isChecked():
						break
					if i < first_fr_num:
						_,_ = cap.read()
					else:
						progress_area.setValue(int((i - first_fr_num) / (last_fr_num - first_fr_num + 1) * 100)) # progressbar
						try:
				#if 1:
							_,test_frame = cap.read()
							test_frame = cv2.cvtColor(test_frame, cv2.COLOR_BGR2GRAY)
							if self.ui.checkBox_Rotate.isChecked():
								test_frame = ndimage.rotate(test_frame, float(self.ui.RotationAngle.text()))
#######
							test_frame = test_frame[crop_y0:crop_y1, crop_x0:crop_x1]
							test_frame = cv2.resize(test_frame, None, fx=scale, fy=scale, interpolation = cv2.INTER_CUBIC)
#######
						#if 1:
							try:
								labels = oneshot_areas_select(test_frame, self.savepath, "Proc",
									borders_remove1, borders_remove2,
									start_coeff2 = start_coeff2,
									min_area2 = min_area2,
									step_coeff2 = step_coeff2,
									sm_blur_val = sm_blur_val,
									start_coeff1 = start_coeff1,
									min_area1 = min_area1,
									step_coeff1 = step_coeff1,
									plt_all = False) #tmp_thr=float(self.ui.lineEdit_Thr2.text())
								data_areas = oneshot_areas_proc(labels, test_frame, self.savepath, plt_all=False)
								print(data_areas)
								df_tmp = pd.DataFrame(data_areas, index=[i])
								if df_data_areas.empty:
									df_data_areas = df_tmp
								else:
									df_data_areas = pd.concat([df_data_areas, df_tmp])

							except:
								print(f"Process failed for area {i}") #Probably we have to create a list
						except:
							print(f'Bad frame {i}')
					i += 1
				cap.release()
				if not self.ui.pushButton_BreakSegm.isChecked():
					progress_area.setValue(100)
				self.ui.pushButton_StartSegm.setChecked(False)
				self.ui.pushButton_BreakSegm.setChecked(False)
				self.ui.pushButton_BreakSegm.setEnabled(False)
				for q in [self.ui.tab1,self.ui.tab3,self.ui.tab4,self.ui.tab_3,self.ui.lineEdit_AreaScale,
						self.ui.lineEdit_14,self.ui.lineEdit_15,
						self.ui.lineEdit_TestFrN,self.ui.lineEdit_sm_blur_val,self.ui.pushButton_10,self.ui.pushButton_11,
						self.ui.lineEdit_minarea_1,self.ui.lineEdit_startval1,self.ui.lineEdit_stepcoeff1,
						self.ui.lineEdit_minarea_2,self.ui.lineEdit_startval2,self.ui.lineEdit_stepcoeff2,
						#self.ui.horizontalSlider,
						self.ui.pushButton_TestCall,self.ui.pushButton_StartSegm,
						self.ui.pushButton_a1_b,self.ui.pushButton_a1_t,self.ui.pushButton_a1_l,self.ui.pushButton_a1_r,
						self.ui.pushButton_a2_b,self.ui.pushButton_a2_t,self.ui.pushButton_a2_l,self.ui.pushButton_a2_r,
						self.ui.pushButton_crop_2, self.ui.lineEdit_30, self.ui.checkBox_Rotate_2, self.ui.RotationAngle_2]:

					q.setEnabled(True)

				# change column order
				df_data_areas = df_data_areas[['substr_angle', 'tip_sample_contact_length', 'dx_tip_sample', 'sample_surface_contact_length', 'dx_sample_surface', 'mid_projection_x', 'mid_projection_y', 'len_mid', 'square_x', 'square_y', 'ellipse_major', 'ellipse_minor', 'area', 'orientation']]
				# df_data_areas stores frame number as index
				df_data_areas.to_csv(f"{self.savepath}_areas.csv", index=True, index_label='frame_number')
				df_data_areas = None
			except:
				QMessageBox.about(self, "Unexpected error", "Please, report to developers")


	def btnClicked_FullContoursAnalysis(self):
		plt.close()


		a = self.ui.lineEdit_18.text()
		b = self.ui.lineEdit_17.text()
		ln = self.ui.lineEdit_12.text()
		angle = self.ui.RotationAngle.text()

		todo = False
		try:
			a, b, ln = int(a), int(b), int(ln)
			if a > b:
				a,b = b,a

			if a < 0: a = 0
			if b < 0: b = 0
			if a > ln: a = ln
			if b > ln: b = ln

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame numbers")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
			todo = False

		try:
			scale = float(self.ui.lineEdit_AreaScale.text())
		except:
			QMessageBox.about(self, "Input error", "Please, check scale")
			todo = False

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo:

			for q in [self.ui.tab1,self.ui.tab3,self.ui.tab_2,self.ui.tab,self.ui.tab4
				,self.ui.pushButton_12,self.ui.pushButton_13,self.ui.checkBox_8
				,self.ui.lineEdit_TestFrN_2,self.ui.lineEdit_AreaScale_2
				,self.ui.pushButton_TestCall_2,self.ui.pushButton_StartSegm_2
				,self.ui.lineEdit_18,self.ui.lineEdit_17,self.ui.pushButton_crop_3
				,self.ui.lineEdit_31, self.ui.checkBox_Rotate_3 ,self.ui.RotationAngle_3]:
				q.setEnabled(False)
			self.ui.pushButton_BreakSegm_2.setEnabled(True)
			self.ui.pushButton_BreakSegm_2.setChecked(False)

			try:

				first_fr_num = a
				last_fr_num = b
				cap = cv2.VideoCapture(self.name_of_videofile)
				len_cap = ln

				area_output = open(self.savepath+'_contours.shpx','w')
				area_output.write('Parameters:\n')
				if self.ui.checkBox_Rotate.isChecked():
					area_output.write('Rotation:\t'+str(self.ui.RotationAngle.text())+'\n')
				area_output.write('Contours scale:\t'+self.ui.lineEdit_AreaScale_2.text()+'\n')
				#area_output.write('Scale:\t'+self.ui.lineEdit_AreaScale_2.text()+'\n')
				area_output.write(self.separator+'\n')
				area_output.write('Frame number'+'\t'
							+'Tip-Sample'+'\t' + 'Sample-Surface'+'\t'
							+'mid_x_proj'+'\t'+ 'mid_y_proj'+'\t'
							+'mid_len' +'\t' + 'Volume_by_top'+'\t' +
							'Volume_by_bottom'+'\t' + 'invvol_t'+'\t' +
							'invvol_b'+'\t' + 'corrvol_top'+'\t' +
							'corrvol_b'+'\n')

				i=0
				while i <= last_fr_num:
					self.ui.lineEdit_CurrentFrameAreas_2.setText(str(i))
					QApplication.processEvents()
					if self.ui.pushButton_BreakSegm_2.isChecked():
						break
					if i < first_fr_num:
						_,_=cap.read()
					else:
						self.ui.progressBar_3.setValue((i-first_fr_num)/(last_fr_num-first_fr_num+1)*100) # progressbar
						try:
							_,test_frame = cap.read()
							test_frame = cv2.cvtColor(test_frame, cv2.COLOR_BGR2GRAY)
#######
							if self.ui.checkBox_Rotate.isChecked():
								test_frame = ndimage.rotate(test_frame, float(self.ui.RotationAngle.text()))
							test_frame = test_frame[crop_y0:crop_y1,crop_x0:crop_x1]
							test_frame=cv2.resize(test_frame, None,fx=scale, fy=scale, interpolation = cv2.INTER_CUBIC)

#					try:
							data_areas = oneshot_contours(test_frame,self.savepath,plt_all=False)
							print(data_areas)
							area_output.write(str(i)+'\t'
								+str(data_areas['ts_len'])+'\t' + str(data_areas['sp_len'])+'\t'
								+str(data_areas['mid_x'])+'\t'+ str(data_areas['mid_y'])+'\t'+
								str(data_areas['mid_len']) +'\t' + str(data_areas['vol_t'])+'\t' +
								str(data_areas['vol_b'])+'\t' + str(data_areas['invvol_t'])+'\t' +
								str(data_areas['invvol_b'])+'\t' + str(data_areas['corrvol_t'])+'\t' +
								str(data_areas['corrvol_b'])+'\n')
						except:
							print('bad frame', i)

					i+=1
				cap.release()
				if not self.ui.pushButton_BreakSegm_2.isChecked():
					self.ui.progressBar_3.setValue(100)

				for q in [self.ui.tab1,self.ui.tab3,self.ui.tab_2,self.ui.tab,self.ui.tab4
						,self.ui.pushButton_12,self.ui.pushButton_13,self.ui.checkBox_8
						,self.ui.lineEdit_TestFrN_2,self.ui.lineEdit_AreaScale_2
						,self.ui.pushButton_TestCall_2,self.ui.pushButton_StartSegm_2
						,self.ui.lineEdit_18,self.ui.lineEdit_17,self.ui.pushButton_crop_3
						,self.ui.lineEdit_31, self.ui.checkBox_Rotate_3 ,self.ui.RotationAngle_3]:
					q.setEnabled(True)
				self.ui.pushButton_BreakSegm_2.setEnabled(False)
				self.ui.pushButton_BreakSegm_2.setChecked(False)
				self.ui.pushButton_StartSegm_2.setChecked(False)

				area_output.close()
			except:
				QMessageBox.about(self, "Unexpected error", "Please, report to developers")

	def btnClicked_TestCall(self):
		plt.close()
		#global allow_contour
		#global allow_area

		a = self.ui.lineEdit_TestFrN.text()
		ln = self.ui.lineEdit_12.text()
		angle = self.ui.RotationAngle.text()

		todo = False
		try:
			a, ln = int(a), int(ln)
			if a < 0: a = 0
			if a > ln: a = ln

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame number")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
			todo = False

		try:
			start_coeff2=float(self.ui.lineEdit_startval2.text())
			min_area2=int(self.ui.lineEdit_minarea_2.text())
			step_coeff2=float(self.ui.lineEdit_stepcoeff2.text())
			sm_blur_val = int(self.ui.lineEdit_sm_blur_val.text())
			start_coeff1=float(self.ui.lineEdit_startval1.text())
			min_area1=int(self.ui.lineEdit_minarea_1.text())
			step_coeff1=float(self.ui.lineEdit_stepcoeff1.text())
			scale = float(self.ui.lineEdit_AreaScale.text())
		except:
			QMessageBox.about(self, "Input error", "Please, check parameters")
			todo = False

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo:

			old_imgs = ['initial', 'dist1', 'marker1', 'rw1', 'init2', 'dist2', 'markers2', 'rw2', 'particle', 'fin']
			for i in old_imgs:
				try:
					os.remove(self.savepath +'_Test_area_'+i+'.png')
				except:
					print('No image to remove')

			#create arrays for two sets of "Preserve borders" buttons
			borders_remove1 = []
			borders_remove2 = []
			for i in ['l','r','b','t']:
				borders_remove1.append(getattr(self.ui,"pushButton_a1_"+i).isChecked())
				borders_remove2.append(getattr(self.ui,"pushButton_a2_"+i).isChecked())
			if 1:
			#try:
				test_frame = get_n_frame(self.name_of_videofile, 0, ln, a)
				test_frame = cv2.cvtColor(test_frame, cv2.COLOR_BGR2GRAY)
#######
				if self.ui.checkBox_Rotate.isChecked():
					test_frame = ndimage.rotate(test_frame, angle)
				print("crop",crop_y0,crop_y1,crop_x0,crop_x1)
				test_frame = test_frame[crop_y0:crop_y1,crop_x0:crop_x1]
				test_frame = cv2.resize(test_frame, None,fx=scale, fy=scale, interpolation = cv2.INTER_CUBIC)
#######
				labels = oneshot_areas_select(test_frame,self.savepath,"Test", borders_remove1, borders_remove2,
								start_coeff2=start_coeff2,
								min_area2=min_area2,
								step_coeff2=step_coeff2,
								sm_blur_val=sm_blur_val,
								start_coeff1=start_coeff1,
								min_area1=min_area1,
								step_coeff1=step_coeff1,
								plt_all=True)
				data_areas = oneshot_areas_proc(labels,test_frame,self.savepath,plt_all=True)
				print(data_areas)
			#except:
			#	QMessageBox.about(self, "Error", "Test failed for areas")

			try:
				a = ['imgIni', 'imgDist1', 'imgMark1', 'imgRw1', 'imgCut', 'imgDist2', 'imgMark2', 'imgRw2', 'imgParticle', 'imgFin']
				b = ['initial', 'dist1', 'marker1', 'rw1', 'init2', 'dist2', 'markers2', 'rw2', 'particle', 'fin']
				i = 0
				while i < len(a):
					getattr(self.ui,'label_'+a[i]).setPixmap(
							QPixmap(self.savepath +'_Test_area_'+b[i]+'.png').scaled(
							getattr(self.ui,'label_'+a[i]).size(),
							QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
					getattr(self.ui,'label_'+a[i]).show()
					i += 1
			except:
					QMessageBox.about(self, "Error", "Error in the visualisation of results")

		self.ui.pushButton_TestCall.setChecked(False)

	def btnClicked_TestCall_contours(self):
		plt.close()

		a = self.ui.lineEdit_TestFrN_2.text()
		ln = self.ui.lineEdit_12.text()
		angle = self.ui.RotationAngle.text()

		todo = False
		try:
			a, ln = int(a), int(ln)
			if a < 0: a = 0
			if a > ln: a = ln

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame number")

		try:
			angle = float(angle)
		except:
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")
			todo = False

		try:
			scale = float(self.ui.lineEdit_AreaScale.text())
		except:
			QMessageBox.about(self, "Input error", "Please, check scale")
			todo = False
		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		#print(self.name_of_videofile, test_fr_num, allow_area, allow_contour)
		if todo:
			try:
				test_frame = get_n_frame(self.name_of_videofile,0,ln, a)
				test_frame = cv2.cvtColor(test_frame, cv2.COLOR_BGR2GRAY)
#######
				if self.ui.checkBox_Rotate.isChecked():
					test_frame = ndimage.rotate(test_frame, angle)
				test_frame = test_frame[crop_y0:crop_y1,crop_x0:crop_x1]
				test_frame=cv2.resize(test_frame, None,fx=scale, fy=scale, interpolation = cv2.INTER_CUBIC)
#######
				try:
					os.remove(self.savepath +'_contours2.png')
				except:
					pass

				try:
					_ = oneshot_contours(test_frame,self.savepath,plt_all=True)
				except:
					QMessageBox.about(self, "Error", "Test failed for contours")
				#'''
				try:
					self.ui.label_imgIni_2.setPixmap(QPixmap(self.savepath +'_c_initial.png').scaled(self.ui.label_imgIni_2.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
					self.ui.label_imgIni_2.show()
					self.ui.label_imgRw2_2.setPixmap(QPixmap(self.savepath +'_contours2.png').scaled(self.ui.label_imgRw2_2.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
					self.ui.label_imgRw2_2.show()
				except:
					QMessageBox.about(self, "Error", "Error in the visualisation of results")
		#'''
				self.ui.pushButton_TestCall.setChecked(False)
			except:
				QMessageBox.about(self, "Error", "Error in the test frame")



	def pnm_recall(self,t):
		global scalebar_length_pix
		global pnm

		sizes = {'nm':1,'μm':10**3,'mm':10**6,'m':10**9}

		if t=="Scale":
			scale = sizes[self.ui.comboBox.currentText()]
			try:
				pnm = float(self.ui.lineEdit_2.text())/scalebar_length_pix*scale
			except:
				pnm = None
		else:
			try:
				mag = float(self.ui.lineEdit.text())
				sc = float(self.ui.lineEdit_8.text())
				pnm = sc/mag
			except:
				pnm = None

	def radiobutton_check_Young(self):
		if self.ui.radioButton_Y_lin.isChecked():
			self.ui.k_tip.setEnabled(True)
			for i in ['probe_modulus','substr_modulus','probe_poiss','substr_poiss']:
				getattr(self.ui,i).setEnabled(False)
			self.ui.tip_radius.setEnabled(False)
			self.ui.cyl_radius.setEnabled(False)
			self.ui.tip_radius_2.setEnabled(False)
			self.ui.sphere_radius.setEnabled(False)
			
		if self.ui.radioButton_Y_Sph.isChecked():
			self.ui.k_tip.setEnabled(False)
			for i in ['probe_modulus','substr_modulus','probe_poiss','substr_poiss']:
				getattr(self.ui,i).setEnabled(True)
			self.ui.tip_radius.setEnabled(False)
			self.ui.cyl_radius.setEnabled(False)
			self.ui.tip_radius_2.setEnabled(True)
			self.ui.sphere_radius.setEnabled(True)
		if self.ui.radioButton_Y_Cyl.isChecked():
			self.ui.k_tip.setEnabled(False)
			for i in ['probe_modulus','substr_modulus','probe_poiss','substr_poiss']:
				getattr(self.ui,i).setEnabled(True)
			self.ui.tip_radius.setEnabled(True)
			self.ui.cyl_radius.setEnabled(True)
			self.ui.tip_radius_2.setEnabled(False)
			self.ui.sphere_radius.setEnabled(False)

	def radiobutton_check(self):

		global scalebar_length_pix

		if self.ui.radioButton_pnmbybar.isChecked():
			self.ui.lineEdit_2.setEnabled(True)
			self.ui.comboBox.setEnabled(True)
			self.ui.lineEdit.setEnabled(False)
			self.ui.lineEdit_8.setEnabled(False)

			self.pnm_recall("Scale")

		if self.ui.radioButton_pnmbymag.isChecked():
			self.ui.lineEdit_2.setEnabled(False)
			self.ui.comboBox.setEnabled(False)
			self.ui.lineEdit.setEnabled(True)
			self.ui.lineEdit_8.setEnabled(True)
			self.pnm_recall("Mag")
		global pnm
		s = str(pnm)
		if len(s)>7:
			s=str(round(pnm,5))
		if pnm:
			self.ui.lineEdit_pnm.setText(s)
			self.ui.lineEdit_pnm2.setText(s)
	'''
	def checkBox_ArC(self):
		global allow_contour
		global allow_area
#		allow_area = self.ui.checkBox_Area.isChecked()
#		allow_contour = self.ui.checkBox_Cont.isChecked()
		if allow_area or allow_contour:
			self.ui.lineEdit_Thr2_2.setEnabled(True)
			self.ui.lineEdit_Thr2.setEnabled(True)
			self.ui.lineEdit_AreaScale.setEnabled(True)
			self.ui.lineEdit_TestFrN.setEnabled(True)
			self.ui.pushButton_TestCall.setEnabled(True)
		else:
			self.ui.lineEdit_Thr2_2.setEnabled(False)
			self.ui.lineEdit_Thr2.setEnabled(False)
			self.ui.lineEdit_AreaScale.setEnabled(False)
			self.ui.lineEdit_TestFrN.setEnabled(False)
			self.ui.pushButton_TestCall.setEnabled(False)
	'''
	def checkBox_allow_out(self):
		if not self.ui.checkBox_3.isChecked():
			self.ui.pushButton_9.setEnabled(True)
		elif self.ui.checkBox_3.isChecked():
			self.ui.pushButton_9.setEnabled(False)

	def btnClicked_import_hysdata(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.hys_filename, _ = QFileDialog.getOpenFileName(self,"Select hysitron data file", os.getcwd(),"txt files (*.hystxt);;txt files (*.txt);;All Files (*)", options=options)
		self.ui.lineEdit_hys_data_import.setText(self.hys_filename)

	def btnClicked_import_olddata(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.old_filename, _ = QFileDialog.getOpenFileName(self,"Select the data file", os.getcwd(),"txt files (*.txt);;txt files (*.txt);;All Files (*)", options=options)
		self.ui.lineEdit_video_data_import_2.setText(self.old_filename)

	def btnClicked_import_videodata(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.video_filename, _ = QFileDialog.getOpenFileName(self,"Select video data file", os.getcwd(), options=options)
		self.ui.lineEdit_video_data_import.setText(self.video_filename)

	def btnClicked_import_corrdata(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.corr_filename, _ = QFileDialog.getOpenFileName(self,"Select video data file", os.getcwd(), options=options)
		self.ui.lineEdit_corr_data_import.setText(self.corr_filename)

	def btnClicked_AreasOpen(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.area_data_import, _ = QFileDialog.getOpenFileName(self,"Select areas data file", os.getcwd(), options=options)
		self.ui.lineEdit_26.setText(self.area_data_import)
		self.ui.pushButton_ConvertAreas.setEnabled(True)

	def btnClicked_ContoursOpen(self):
		plt.close()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		self.contours_data_import, _ = QFileDialog.getOpenFileName(self,"Select contours data file", os.getcwd(), options=options)
		self.ui.lineEdit_27.setText(self.contours_data_import)
		self.ui.pushButton_ConvertContours.setEnabled(True)


	def update_crop(self,fr):
		#fr = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)
		global init_img_width, init_img_height
		global crop_x0,crop_x1,crop_y0,crop_y1
		init_img_width = len(fr[0])
		init_img_height = len(fr)
		crop_x0,crop_x1,crop_y0,crop_y1 = 0,init_img_width,0,init_img_height
		#print(crop_x0,crop_x1,crop_y0,crop_y1)

	def btnClicked_selectfile(self):
		plt.close()
		self.areas_recall()
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getOpenFileName(self,"Select a video file", os.getcwd(),"Video files (*.avi);;All Files (*)", options=options)

		############## Is the file valid?
		###############
		if fileName != '':
			self.name_of_videofile = fileName
			print(self.name_of_videofile)

			try:
				#Probe the img
				cap = cv2.VideoCapture(self.name_of_videofile)
				len_cap = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
				print('Video lengh',len_cap)
				i = 0
				while i < len_cap:
					try:
						_,img = cap.read()
						break
					except:
						print(i)
					i+=1
				cap.release()

				self.update_crop(img)

				self.ui.lineEdit_3.setText(self.name_of_videofile)
				#self.ui.lineEdit_6.setText(str(int(cv2.VideoCapture(self.ui.lineEdit_3.text()).get(cv2.CAP_PROP_FRAME_COUNT))))
				self.ui.checkBox_3.setEnabled(True)
				self.ui.checkBox_7.setEnabled(True)
				self.ui.checkBox_8.setEnabled(True)

				for q in ['crop','crop_2','crop_3','SelectScalebar','8','5',
						'TestCall','TestCall_2']:
					getattr(self.ui,"pushButton_"+q).setEnabled(True)

				for q in ['sm_blur_val','14','15','18','17',
						'startval2','minarea_2','stepcoeff2','startval1',
						'minarea_1','stepcoeff1','AreaScale','AreaScale_2',
						'TestFrN','TestFrN_2']:
					getattr(self.ui,"lineEdit_"+q).setEnabled(True)

				if self.ui.checkBox_3.isChecked() is False or self.ui.lineEdit_7.text() == '':
					self.savedir = self.name_of_videofile[:self.name_of_videofile.rfind('/')]
					self.ui.lineEdit_7.setText(self.savedir)

				for q in ['6','5','15','13','17',
						'19','25','16','12']:
					getattr(self.ui,"lineEdit_"+q).setText(str(len_cap))
				for q in ['4','14','18','11']:
					getattr(self.ui,"lineEdit_"+q).setText('0')

				self.ui.checkBox_Rotate_3.setEnabled(True)
				self.ui.checkBox_Rotate_2.setEnabled(True)
				self.ui.checkBox_Rotate.setEnabled(True)
				self.ui.lineEdit_10.setText(str(min(10,len_cap)))
				self.radiobutton_check()
				self.ui.pushButton_StartSegm.setEnabled(True)
				self.ui.pushButton_StartSegm_2.setEnabled(True)

				#video_fullname = self.name_of_videofile
				self.savepath = os.path.join(self.savedir, self.name_of_videofile[self.name_of_videofile.rfind('/')+1 : self.name_of_videofile.rfind('.')])

				for q in ['10','11','9','4','5','20','14','15','18','17',
						'TestFrN','TestFrN_2']:
					getattr(self.ui,"lineEdit_"+q).setValidator(QIntValidator(0,len_cap))

				#global pnm
				#pnm = None
				self.ui.lineEdit_pnm2.setText('1')
			except:
				QMessageBox.about(self, "Error", "Wrong file selected")
		else:
			QMessageBox.about(self, "Error", "Please, select a file")

	#Click on the SaveDir button
	def btnClicked_select_saveDir(self):
		plt.close()
		# select directory to self.save things
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog

		try:
			self.savedir = str(QFileDialog.getExistingDirectory(self, "Select Directory", os.getcwd(), options=options))
			self.ui.lineEdit_7.setText(self.savedir)

			for q in ['3','7','8']:
				getattr(self.ui, f'checkBox_{q}').setEnabled(True)
		except:
			QMessageBox.about(self, "Alarm!", "Wrong dir name")

		#test write access for output dir
		try:
			name_of_testfile = 'mech-video-test-for-write-access'
			path_to_testfile = os.path.join(self.savedir, name_of_testfile)
			if name_of_testfile not in os.listdir(self.savedir):
				with open(path_to_testfile, 'w') as file:
					file.write('test')
				os.remove(path_to_testfile)
				print('Write access OK')
		except:
			QMessageBox.about(self, "Alarm!", "No write access for chosen directory")

	# frames preview


	#If Slider updates in the frame preview window


	# Start DIC
	def btnClicked_start_tracker(self):
		if int(self.ui.lineEdit_6.text()) <= 1:
			QMessageBox.about(self, "Something is wrong with videofile or codec")

		a = self.ui.lineEdit_4.text()
		b = self.ui.lineEdit_5.text()
		c = self.ui.lineEdit_9.text()
		ln = self.ui.lineEdit_6.text()
		angle = self.ui.RotationAngle.text()
		todo = False
		try:
			a, b, c, ln = int(a), int(b), int(c), int(ln)
			if a > b:
				a,b = b,a

			if a < 0: a = 0
			if b < 0: b = 0
			if a > ln: a = ln
			if b > ln: b = ln

			todo = True
		except:
			QMessageBox.about(self, "Input error", "Please, check frame numbers")

		try:
			angle = float(angle)
		except:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check the rotation angle")

		if ln <= 1 and todo:
			todo = False
			QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

		if todo:
			try:# if mag and scale were chosen
				global pnm
				#avg = 1 #Number of frames to be averaged
				#should not be changed yet!

				global start
				global stop
				start, stop, corr_n = a, b, c

##########################################
				##Probably, it is not the correct decision
				#have to be fixed
				corr_n = max(start, corr_n)
				corr_n = min(stop, corr_n)
##########################################

				frame_numbers, all_X, all_Y = tracker(self.name_of_videofile,rectangle_coordinates,
									(crop_x0,crop_x1,crop_y0,crop_y1),
									(start,stop,corr_n),
									rotate = self.ui.checkBox_Rotate.isChecked(),
									rotation_angle = float(self.ui.RotationAngle.text()))
				all_X = np.array(all_X)
				all_Y = np.array(all_Y)
				displ = np.sqrt(all_X ** 2 + all_Y ** 2)

				#if len(output) > 1:
				#	substrate = output[1]
				i = 0
				while i < len(all_X):
					plt.plot(frame_numbers, all_X[i,:]-all_X[i][0], label= 'X')
					plt.plot(frame_numbers, all_Y[i,:]-all_Y[i][0], label= 'Y')
					plt.xlabel("Frame number")
					plt.ylabel("$d$, $pix$")
					plt.legend()
					plt.title("Movement of area #"+str(i+1))
					#plt.self.savefig(self.ui.lineEdit_3.text()[:self.ui.lineEdit_3.text().rfind('.')]+'_mv.png')

					plt.savefig(self.savepath +'_mv'+str(i)+'.png')
					plt.close()
					i += 1

				self.ui.label_10.setPixmap(QPixmap(self.savepath +'_mv0.png').scaled(self.ui.label_10.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
				self.ui.label_10.show()

				self.ui.label_13.setPixmap(QPixmap(self.savepath +'_mv1.png').scaled(self.ui.label_13.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
				self.ui.label_13.show()

				#write output for all areas in pix
				o_full_pix = open(f'{self.savepath}_mv_pix.txt','w')
				if self.ui.checkBox_Rotate.isChecked():
					o_full_pix.write(f'Rotation:\t{self.ui.RotationAngle.text()}\n')
				o_full_pix.write(f'pnm\t{pnm}\n')
				o_full_pix.write(f'{self.separator}\n')
				header_str ="Frame\t"
				i = 1
				while i <= len(all_X):
					header_str += f"Area {i} x , pix\tArea {i} x , pix\tArea {i} displacement, pix\t"
					i += 1
				o_full_pix.write(header_str + '\n')

				for i in range(len(frame_numbers)):
					data_str = str(frame_numbers[i]) + '\t'
					j = 0
					while j < len(all_X):
						#print(all_X[0],displ[0])
						data_str += f"{all_X[j,i]}\t{all_Y[j,i]}\t{displ[j,i]}\t"
						j += 1
					o_full_pix.write(data_str + '\n')
				o_full_pix.close()

				#Conversion to nm
				all_X = all_X * pnm
				all_Y = all_Y * pnm
				displ = displ * pnm

				#write output for all areas in nm
				o_full_nm = open(self.savepath + '_mv_nm.txt','w')
				if self.ui.checkBox_Rotate.isChecked():
					o_full_nm.write(f'Rotation:\t{self.ui.RotationAngle.text()}\n')
				o_full_nm.write(f'pnm\t{pnm}\n')
				o_full_nm.write(self.separator + '\n')
				header_str ="Frame\t"
				i = 1
				while i <= len(all_X):
					header_str += f"Area {i} x , nm\tArea {i} y , nm\tArea {i} displacement, nm\t"
					i += 1
				o_full_nm.write(header_str + '\n')

				for i in range(len(frame_numbers)):
					data_str = str(frame_numbers[i]) + '\t'
					j = 0
					while j < len(all_X):
						data_str += f"{all_X[j,i]}\t{all_Y[j,i]}\t{displ[j,i]}\t"
						j += 1
					o_full_nm.write(data_str + '\n')
				o_full_nm.close()

				progress.setValue(100)
				self.ui.pushButton_start_tracker.setChecked(False)
			except:
				QMessageBox.about(self, "Tracker error", "Please, report the case to developers")

	def btnClicked_P3(self):
		### parse data
		# video
		self.time_v, self.disp_v, self.pnmi = parse_video_data(self.video_filename, self.separator, self.ui.corr_bckgr.isChecked())
		self.time_v, self.disp_v = self.time_v -self.time_v[0], self.disp_v - self.disp_v[0] #Only relative coordinates matter
		# Hys
		self.time_h, self.disp_h, self.load_h = parse_hys_data(self.hys_filename)

		arguments = {}
		arguments['time_v'] = self.time_v
		arguments['disp_v'] = self.disp_v
		arguments['pnmi'] = self.pnmi
		arguments['time_h'] = self.time_h
		arguments['disp_h'] = self.disp_h
		arguments['load_h'] = self.load_h
		
		self.window_P3 = WindowComparison(arguments)
		self.window_P3.run()
		# does not wait for plt window to close! parameters are updated only when btnClicked_P7 is called!
		self.ui.pushButton_7.setEnabled(True)

	def plot_old(self):
	
		old_time_h, old_disp, old_load_h = parse_old_data(self.old_filename)
	
		###################Fig1##################
		fig = plt.figure(figsize= (8, 6))
		grid = plt.GridSpec(20, 4, hspace=0.2, wspace=0.2)
		ax0=plt.subplot(grid[0:8, 0:])
		ax1=plt.subplot(grid[11:20, 0:])
		ax2=ax1.twinx()
		ax0.set_xlabel("$h$, $nm$")
		ax0.set_ylabel("$F$, $\mu N$", color= 'b')
		ax1.set_xlabel("$t$, $s$")
		ax2.set_ylabel("$h$, $nm$", color= 'r')
		ax1.set_ylabel("$F$, $\mu N$", color= 'b')
		ax0.grid(True)
		ax1.grid(True)
		#ax2.grid(True)
		plt.subplots_adjust(right=0.85, left=0.15, top=0.9, bottom=0.13)

		#ax1.axes.get_xaxis().set_visible(False)
		#lim0=1.1*min(min(self.load_h), min(disp))
		#lim1=1.1*max(max(self.load_h), max(disp))
#		k=1.1*max(max(fit_y), max(fit_x))

		for t in [ax0,ax1,ax2]:
			t.xaxis.label.set_size(16)
			t.yaxis.label.set_size(16)
			t.tick_params(labelsize=14)

		ax0.plot(old_disp, old_load_h, "b", label= 'F(d)')

		#legend = ax[0].legend(loc= 'upper left', fontsize=16)
		ax1.plot(old_time_h, old_load_h,"b")
		ax2.plot(old_time_h, old_disp, "r")
		#ax0.set_title(fname)

		#plt.show()
		#plt.savefig('f1_left.png', dpi=600)
		fig.subplots_adjust(hspace=0.3)
		plt.savefig(self.old_filename + '.png',dpi = 1200)
		plt.close(fig)
		plt.close()


	def btnClicked_P7(self):
		self.time_h, self.disp_h, self.time_v, self.disp_v, self.parameters = self.window_P3.get_parameters()
		print(self.parameters)
		_, disp, _, _ = Convert(self.time_h, self.disp_h, self.time_v, self.disp_v, self.parameters)

		v_fullname = self.ui.lineEdit_video_data_import.text()
		if self.savedir != None:
			savedir = self.savedir
		else:
			savedir = v_fullname[:v_fullname.rfind('/')]
		fname = v_fullname[v_fullname.rfind('/')+1:]

		savepath = os.path.join(savedir, v_fullname[v_fullname.rfind('/')+1 : v_fullname.rfind('.')])

		file = open(savepath + '_corrected.txt','w',encoding='utf-8')
#		for row in self.info:
#			file.write(row+'\n')

		###############Wrong value!
		if self.ui.checkBox_Rotate.isChecked():
			file.write('Rotation:\t' + str(self.ui.RotationAngle.text()) + '\n')
		##############
		file.write('Time scale\t' + str(self.parameters['time_scale']) + '\n')
		file.write('Time shift\t' + str(self.parameters['time_shift']) + '\n')
		file.write('Displacement scale\t' + str(self.parameters['disp_scale']) + '\n')
		file.write('Displacement shift\t' + str(self.parameters['disp_shift']) + '\n')
		file.write('a\t' + str(self.parameters['line_scale']) + '\n')
		file.write('sin ampl\t' + str(self.parameters['sine_amplitude']) + '\n')
		file.write('sin freq\t' + str(self.parameters['sine_frequency']) + '\n')
		file.write('sin phase shift\t' + str(self.parameters['sine_shift']) + '\n')
		file.write('pnm\t' + self.pnmi + '\n')
		file.write(self.separator + '\n')
		file.write('Time, s\tDisplacement, nm\tLoad, µN\n')
		print('1')
		for i in range(len(self.time_h)):
			file.write(f"{self.time_h[i]}\t{disp[i]}\t{self.load_h[i]}\n")
		file.close()

		###################Fig1##################
		fig = plt.figure(figsize= (8, 6))
		grid = plt.GridSpec(20, 4, hspace=0.2, wspace=0.2)
		ax0=plt.subplot(grid[0:8, 0:])
		ax1=plt.subplot(grid[11:20, 0:])
		ax2=ax1.twinx()
		ax0.set_xlabel("$h$, $nm$")
		ax0.set_ylabel("$F$, $\mu N$", color= 'b')
		ax1.set_xlabel("$t$, $s$")
		ax2.set_ylabel("$h$, $nm$", color= 'r')
		ax1.set_ylabel("$F$, $\mu N$", color= 'b')
		ax0.grid(True)
		ax1.grid(True)
		#ax2.grid(True)
		plt.subplots_adjust(right=0.85, left=0.15, top=0.9, bottom=0.13)

		#ax1.axes.get_xaxis().set_visible(False)
		#lim0=1.1*min(min(self.load_h), min(disp))
		#lim1=1.1*max(max(self.load_h), max(disp))
#		k=1.1*max(max(fit_y), max(fit_x))

		for t in [ax0,ax1,ax2]:
			t.xaxis.label.set_size(16)
			t.yaxis.label.set_size(16)
			t.tick_params(labelsize=14)

		ax0.plot(disp,self.load_h, "b", label= 'F(d)')

		#legend = ax[0].legend(loc= 'upper left', fontsize=16)
		ax1.plot(self.time_h,self.load_h,"b")
		ax2.plot(self.time_h,disp,"r")
		ax0.set_title(fname)

		#plt.show()
		#plt.savefig('f1_left.png', dpi=600)
		fig.subplots_adjust(hspace=0.3)
		plt.savefig(savepath+'_corrected.png')
		plt.close(fig)
		plt.close()

		self.ui.label_37.setPixmap(QPixmap(savepath +'_corrected.png').scaled(self.ui.label_37.size(), QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation))
		self.ui.label_37.show()
		self.ui.lineEdit_fps.setText(str(self.parameters['time_scale']))
		self.ui.lineEdit_shift.setText(str(self.parameters['time_shift']))

	def btnClicked_openComparison(self):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		t, _ = QFileDialog.getOpenFileName(self,"Select comparison file", os.getcwd(), options=options)
		try:
			d = open(t,'r',encoding='utf-8')
			while 1:
				s = d.readline()
				if s.startswith(self.separator):
					break
				elif s.startswith('Time scale'):
					fps = s.split('\t')
					self.ui.lineEdit_fps.setText(fps[1])
				elif s.startswith('Time shift'):
					sh = s.split('\t')
					self.ui.lineEdit_shift.setText(sh[1])
				elif s.startswith('pnm'):
					ipnm = s.split('\t')
					self.ui.lineEdit_pnm2.setText(ipnm[1])
			data = d.readline()
			data = d.readlines()
			d.close()
			data_v = [s.split('\t') for s in data]
			data = np.array(data_v, dtype= 'float')
			self.fin_ht = data[:,0]
			self.fin_hx = data[:,1]
			self.fin_hf = data[:,2]
			self.ui.lineEdit_29.setText(t)
			#if self.ui.pushButton_ConvertAreas.isChecked() and self.ui.pushButton_ConvertContours.isChecked():
			#	self.ui.pushButton_FinPlot.setEnabled(True)
		except:
			QMessageBox.about(self, "Alarm!", "Something went wrong. Please, try again")

###############
##############
##############
#################
	def btnClicked_P4(self):
		plt.close()
		todo = True		
		v1 = self.ui.substr_poiss.text()
		v2 = self.ui.probe_poiss.text()

		try:
			v1, v2 = float(v1), float(v2)
			if v1 == 1. or v2 == 1.:
				print('err')
				raise ValueError
		except:
			todo = False
			QMessageBox.about(self, "Error", "Please, check values of Poisson's ratio")			

		if todo:
			#if 1:
			try:
				self.time_c, self.disp_c, self.load_c = parse_corr_data(self.corr_filename,self.separator)
				self.load_c = np.array(self.load_c) * float(self.ui.force_scale.text())
				print("Force scale factor is", self.ui.force_scale.text())
				self.load_c = self.load_c.tolist()
				self.determ_borders()
			except:
				QMessageBox.about(self, "Unexpected error", "Please, report to developers")

	def e_determ(self,xf,ff):
		if self.ui.radioButton_Y_lin.isChecked():
			popt,pconv = scipy.optimize.curve_fit(linear,xf,ff)
			perr = np.sqrt(np.diag(pconv))
			print(popt,perr)
			ttt=popt[0]
			E_tip = float(self.ui.k_tip.text()) #45. #N/m
			self.fit_func = lambda x: x/1000. - ttt#mkN
		elif self.ui.radioButton_Y_Sph.isChecked():
			r_sph=float(self.ui.sphere_radius.text())#300.#/10**9#nm
			r_tip=float(self.ui.tip_radius_2.text())#1000.#/10**9
			vSi=float(self.ui.substr_poiss.text())#.22
			E_Si=float(self.ui.substr_modulus.text())/1000.#.165#TPa
			Es_Si=E_Si/(1-vSi**2)
			vD=float(self.ui.probe_poiss.text())#.07
			E_D=float(self.ui.probe_modulus.text())/1000.#1.140#TPa
			Es_D=E_D/(1-vD**2)

			popt,pconv = scipy.optimize.curve_fit(sneddon,xf,ff,[xf[0],ff[0],(ff[-1]-ff[0])/(xf[-1]-xf[0])])
			perr = np.sqrt(np.diag(pconv))
			print(popt[2],perr)
			ttt=popt[2]/(4./3.)#**(2/3)

			Rts=1/(1/r_tip + 1/r_sph)
			Rsp=r_sph
			self.fit_func = lambda x: (1/(Rts*eab(Es_D,x)**2)**(1/3) + 1/(Rsp*eab(Es_Si,x)**2)**(1/3))**(-3/2)-ttt

		elif self.ui.radioButton_Y_Cyl.isChecked():
			r_tube= float(self.ui.cyl_radius.text())#10.#/10**9#nm
			r_tip = float(self.ui.tip_radius.text()) #250.

			vSi=float(self.ui.substr_poiss.text())#.22
			E_Si=float(self.ui.substr_modulus.text())/1000.#.165#TPa
			Es_Si=E_Si/(1-vSi**2)
			vD=float(self.ui.probe_poiss.text())#.07
			E_D=float(self.ui.probe_modulus.text())/1000.#1.140#TPa
			Es_D=E_D/(1-vD**2)

			#E_star = eab(Es_D,x)
			R_star_E = r_tube*r_tip/(2*r_tube + r_tip)
			eps = (r_tip + r_tube)/r_tube
			gamma = 2*(1 + 0.4*np.log(eps))**3/3/(eps + 0.6)/np.power(eps,0.27)
			print("gamma = ",gamma)
			#gamma = 0.5
			#####Simplification for l_tube >> delta_1
			popt,pconv = scipy.optimize.curve_fit(sneddon,xf,ff,[xf[0],ff[0],(ff[-1]-ff[0])/(xf[-1]-xf[0])])
			perr = np.sqrt(np.diag(pconv))
			print(popt[2],perr)
			ttt=popt[2]/(4./3.)/np.sqrt(R_star_E/gamma)#**(2/3)
			self.fit_func = lambda x: eab(x,Es_D) - ttt

		e_star = fsolve(self.fit_func, .05)
		print("E_star ",e_star[0]*1000,"GPa")
		if self.ui.radioButton_Y_lin.isChecked():
			print("Calib coeff:",E_tip/e_star*float(self.ui.force_scale.text()))
		return e_star,0,popt
################
#################
	def btnClicked_show(self):
		plt.close()

		if self.ui.lineEdit_3.text():

			a = self.ui.lineEdit_9.text()
			ln = self.ui.lineEdit_6.text()
			angle = self.ui.RotationAngle.text()
			todo = False
			try:
				a, ln = int(a), int(ln)
				if a < 0: a = 0
				if a > ln: a = ln
				todo = True
			except:
				QMessageBox.about(self, "Input error", "Please, check frame number")
			try:
				angle = float(angle)
			except:
				todo = False
				QMessageBox.about(self, "Input error", "Please, check the rotation angle")

			if ln <= 1 and todo:
				todo = False
				QMessageBox.about(self, "Input error", "Please, check video length and codec in use")

			if todo:
				global crop_area_select
				crop_area_select = False
				if self.ui.checkBox_select_crop.isChecked():
					crop_area_select = True

				try:
					global img
					img = get_n_frame(self.name_of_videofile,0,ln,a)
					if self.ui.checkBox_Rotate.isChecked():
						img = ndimage.rotate(img, angle)
					if crop_area_select:
						img=img[crop_y0:crop_y1,crop_x0:crop_x1]
					img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)


#####To move as a function of successfully selected areas
					self.ui.pushButton_start_tracker.setEnabled(True)
#####

					global fig
					global current_ax
					fig, current_ax = plt.subplots()

					plt.imshow(img.copy(), cmap='gray')

					global rectangle_coordinates
					rectangle_coordinates = []

					toggle_selector.RS = RectangleSelector(current_ax, line_select_callback_areas,
								   drawtype= 'box', useblit=True,
								   button= [1, 3],  # don't use middle button
								   minspanx=5, minspany=5,
								   spancoords= 'pixels',
								   interactive=True)
					current_ax.set_title('Select the tip and press "Enter"')
					plt.connect('key_press_event',toggle_selector)
					current_ax.axis('off')
					plt.show()
				except:
					QMessageBox.about(self, "Error", "Some problem detected with the selected frame.\nPlease, select another one")
		else:
			QMessageBox.about(self, "Alarm!", "Select a video file!")


	def determ_borders(self):
		#rewrite to new variables
		global determ_d, determ_force, determ_time
		determ_time=self.time_c
		determ_d=self.disp_c
		determ_force=self.load_c
		#plot
		global fig_young
		fig_young = plt.figure(figsize= (8, 6))
		grid = plt.GridSpec(20, 4, hspace=0.2, wspace=0.2)
		ax1=plt.subplot(grid[0:8, 0:])

		#draw Force vs Displacement
		#ax1=plt.subplot(2,1,1)
		global dfplot, fitplot
		dfplot, =plt.plot(determ_d, determ_force, lw=2, color= 'blue')
		fitplot, =plt.plot(determ_d, determ_force, lw=0.25, color= 'red')
		global lp1, rp1
		lp1, = plt.plot(determ_d[1], determ_force[1], "ro")
		rp1, = plt.plot(determ_d[2], determ_force[2], "go")

		#draw time dependency
		#ax2=plt.subplot(2,1,2)
		ax2=plt.subplot(grid[11:20, 0:])
		plt.plot(determ_time, determ_d, lw=2, color= 'blue')
		global lp2, rp2
		lp2, = plt.plot(determ_time[1], determ_d[1], "ro")
		rp2, = plt.plot(determ_time[2], determ_d[2], "go")



		ax1.set_xlabel("$d$, $nm$")
		ax1.set_ylabel("$F$, $\mu N$")
		ax2.set_xlabel("$t$, $s$")
		ax2.set_ylabel("$d$, $nm$")

		plt.subplots_adjust(left=0.1, bottom=0.3,top=.98)
		#draw buttons and sliders
		global left_border_sl, right_border_sl
		left_border = plt.axes([0.15, 0.15, 0.7, 0.025])
		right_border = plt.axes([0.15, 0.1, 0.7, 0.025])
		try:
			left_border_sl = Slider(left_border, 'Left', 1, len(determ_time)-1, valinit=1, valstep=1)
			right_border_sl = Slider(right_border, 'Right', 1, len(determ_time)-1, valinit=2, valstep=1)
		except:
			left_border_sl = Slider(left_border, 'Left', 1, len(determ_time)-1, valinit=1)
			right_border_sl = Slider(right_border, 'Right', 1, len(determ_time)-1, valinit=2)

		#borders update
		left_border_sl.on_changed(update_determ)
		right_border_sl.on_changed(update_determ)
		#recall
		resetax = plt.axes([0.75, 0.025, 0.1, 0.04])
		global Estar_text
		Estar_text = ax1.text(.95,.95,'', size=18, color='red', horizontalalignment='right',
						verticalalignment='top', transform=ax1.transAxes)
		global reset_but
		reset_but = Button(resetax, 'Reset', hovercolor= '0.975')
		reset_but.on_clicked(reset_determ)
		#Young calc
		youngax = plt.axes([0.2, 0.025, 0.1, 0.04])
		global young_button
#       	yu_modulus=0
		young_button = Button(youngax, 'Calc', hovercolor= '0.975')
		young_button.on_clicked(self.young_determ)
		plt.show()

		s1 = int(left_border_sl.val)
		s2 = int(right_border_sl.val)

		return determ_time[s1:s2],determ_d[s1:s2],determ_force[s1:s2]

	def young_determ(self,val):
		Estar_text.set_text('')
		s1 = int(left_border_sl.val)
		s2 = int(right_border_sl.val)
		lp1.set_data(determ_d[s1], determ_force[s1])
		rp1.set_data(determ_d[s2], determ_force[s2])
		lp2.set_data(determ_time[s1], determ_d[s1])
		rp2.set_data(determ_time[s2], determ_d[s2])

		estar,_,popt=self.e_determ(determ_d[s1:s2], determ_force[s1:s2])
		#print(popt)
		
		dfplot.set_data(determ_d[s1:s2], determ_force[s1:s2])
		if self.ui.radioButton_Y_Sph.isChecked() or self.ui.radioButton_Y_Cyl.isChecked():
			fitplot.set_data(determ_d[s1:s2], sneddon(determ_d[s1:s2],popt[0],popt[1],popt[2]))
		if self.ui.radioButton_Y_lin.isChecked():
			fitplot.set_data(determ_d[s1:s2], linear(np.array(determ_d[s1:s2]),popt[0],popt[1]))
		Estar_text.set_text('$E^{*} = $'+str(float(int(estar[0]*10000))/10)+' $GPa$')
		fig_young.canvas.draw_idle()

	def closeEvent(self, event):
		self.ui.pushButton_BreakSegm.setChecked(True)
		self.ui.pushButton_BreakSegm_2.setChecked(True)
		plt.close()
		print('Close')
		event.accept()
	###


def tracker(name_of_videofile, coords_list, crop, frames_n, rotate=False, rotation_angle = 0): #last frames are empty and cause error
	cap = cv2.VideoCapture(name_of_videofile)
	print(len(coords_list),coords_list)
	frame_numbers = []
	border_array = []
	all_X = []
	all_Y = []
	for q in coords_list:
		border_array.append([int(x) for x in q])
		all_X.append([])
		all_Y.append([])
	start, stop, initial = frames_n
	crop_x0,crop_x1,crop_y0,crop_y1 = crop

	frame = get_n_frame(name_of_videofile,start, stop, initial)
	tr_img = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
	#tr_img = tr_img[crop_y0:crop_y1,crop_x0:crop_x1] #Areas are selected from uncutted image
	if rotate:
		tr_img = ndimage.rotate(tr_img, rotation_angle)
	crop_img = []
	for k in border_array:
		crop_img.append(tr_img[k[1]:k[3], k[0]:k[2]])
	print("Areas for DIC were collected")
	#Probably, this part has to be moved to the selection module

	cap = cv2.VideoCapture(name_of_videofile)

	i=0
	while i < stop:
		if i < start:
			_, _ = cap.read()
		else:
			try:
				progress.setValue(int((i-start)/(stop-start)*100)) # progressbar
				_,frame = cap.read()
				frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
				if rotate:
					frame = ndimage.rotate(frame, rotation_angle)
				cr_test = frame[crop_y0:crop_y1,crop_x0:crop_x1]
				#'''
				j = 0
				while j < len(border_array):
					corr_cv = cv2.matchTemplate(cr_test, crop_img[j], cv2.TM_SQDIFF_NORMED) #image,templ,method,result
					_, _, ml, _ = cv2.minMaxLoc(corr_cv) #min_val,max_val,min_loc,max_loc
					#retval, response = cv2.phaseCorrelate(cr_test, crop_img)
					#print(retval,response)
					all_X[j].append(ml[0])
					all_Y[j].append(ml[1])
					j += 1
				frame_numbers.append(i)
			except:
				print('bad frame', i)
				#break
		i += 1
	cap.release()
	return frame_numbers, all_X, all_Y

######Few functions for interactive matplotlib graphs
######
def update_determ(val):
	s1 = int(left_border_sl.val)
	s2 = int(right_border_sl.val)
	lp1.set_data(determ_d[s1], determ_force[s1])
	rp1.set_data(determ_d[s2], determ_force[s2])
	lp2.set_data(determ_time[s1], determ_d[s1])
	rp2.set_data(determ_time[s2], determ_d[s2])
	dfplot.set_data(determ_d[s1:s2], determ_force[s1:s2])
	fitplot.set_data(determ_d[s1:s2], determ_force[s1:s2])
	fig_young.canvas.draw_idle()

def reset_determ(event):
	left_border_sl.reset()
	right_border_sl.reset()

def toggle_selector(event):
	if event.key in ['q']:
		toggle_selector.RS.set_active(False)
		plt.close()

	if event.key in ['enter']:
		#print(x1,x2,y1,y2)
		toggle_selector.RS.set_active(False)
		toggle_selector.RS.set_active(True)
		#print(rectangle_coordinates)
		#print(' RectangleSelector deactivated.')
		if not rectangle_coordinates:
			if crop_area_select:
				rectangle_coordinates.append([ar_x1+crop_x0,ar_y1+crop_y0,ar_x2+crop_x0,ar_y2+crop_y0])
			else:
				rectangle_coordinates.append([ar_x1, ar_y1, ar_x2, ar_y2])
			current_ax.set_title('Now select the substrate and press "Enter"')
			print("Draw another rectangle or press q to exit. ")
		else:
			if crop_area_select:
				rectangle_coordinates.append([ar_x1+crop_x0,ar_y1+crop_y0,ar_x2+crop_x0,ar_y2+crop_y0])
			else:
				rectangle_coordinates.append([ar_x1, ar_y1, ar_x2, ar_y2])
			current_ax.set_title('Select additional area and press "Enter" or press q to exit.')
	if event.key in ['r']:
		toggle_selector.RS.set_active(False)
		toggle_selector.RS.set_active(True)
		global radius
		radius = [abs(ar_x1-ar_x2)/2, abs(ar_y1-y2)/2]
		print('Radius: ',radius)
	plt.draw()



def line_select_callback_crop(eclick, erelease):
	'First function for drawing a rectangle; eclick and erelease are the press and release events'
	global cr_x1, cr_x2, cr_y1, cr_y2
	cr_x1, cr_y1 = eclick.xdata, eclick.ydata
	cr_x2, cr_y2 = erelease.xdata, erelease.ydata

def line_select_callback_scbar(eclick, erelease):
	'First function for drawing a rectangle; eclick and erelease are the press and release events'
	global sc_x1, sc_x2, sc_y1, sc_y2
	sc_x1, sc_y1 = eclick.xdata, eclick.ydata
	sc_x2, sc_y2 = erelease.xdata, erelease.ydata

def line_select_callback_areas(eclick, erelease):
	'First function for drawing a rectangle; eclick and erelease are the press and release events'
	global ar_x1, ar_x2, ar_y1, ar_y2
	ar_x1, ar_y1 = eclick.xdata, eclick.ydata
	ar_x2, ar_y2 = erelease.xdata, erelease.ydata


def toggle_selector_crop(event):
	global crop_x0,crop_x1,crop_y0,crop_y1
	global cr_x1, cr_x2, cr_y1, cr_y2
	if event.key in ['q']:
		toggle_selector_crop.RS.set_active(False)
		plt.close()
	if event.key in ['enter']:
		#print(cr_x1,cr_x2,cr_y1,cr_y2)
		toggle_selector_crop.RS.set_active(False)
		toggle_selector_crop.RS.set_active(True)
		crop_x0,crop_x1,crop_y0,crop_y1 = max(0,int(cr_x1)),min(init_img_width,int(cr_x2)),max(0,int(cr_y1)),min(init_img_height,int(cr_y2))
		print("Crop area changed",crop_x0,crop_x1,crop_y0,crop_y1)
		current_ax3.set_title('Reselect the crop area and press "Enter" or press q to exit.')
	plt.draw()

def toggle_selector_bar(event):
	global scalebar_length_pix
	if event.key in ['q']:
		toggle_selector_bar.RS.set_active(False)
		plt.close()
	if event.key in ['enter']:
		while len(current_ax_scalebar.lines) >= 1:
			del current_ax_scalebar.lines[0]
		#del selected_contour
		toggle_selector_bar.RS.set_active(False)
		toggle_selector_bar.RS.set_active(True)
		x = [int(z) for z in [sc_x1, sc_y1, sc_x2, sc_y2]]
		contours = measure.find_contours(img[x[1]:x[3], x[0]:x[2]], 10)

		scalebar_length_pix = 0
		for n, contour in enumerate(contours):
			width = max(contour[:, 1]) - min(contour[:, 1])
			if width > scalebar_length_pix:
				scalebar_length_pix = width
				scale_n = n
				#print(scalebar_length_pix)
		contour = contours[scale_n]
		selected_contour, = current_ax_scalebar.plot(np.array(contour[:, 1]) + x[0],
							np.array(contour[:, 0]) + x[1], linewidth=2, color='r')
		print("Scalebar length, pix = ", scalebar_length_pix)
		#scalebar_length_pix = max(contour[:, 1]) - min(contour[:, 1])
		current_ax_scalebar.set_title('Reselect the scale bar and press "Enter" or press q to exit.')
	plt.draw()
