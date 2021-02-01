#!/usr/bin/python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import cv2

def frames_preview(f,start,stop,len_cap,r,angle,cr,cr_val):
	plt.close()

	global st
	st = start

	crop_x0, crop_x1, crop_y0, crop_y1 = cr_val
	cap = cv2.VideoCapture(f)
	images= []
	i=0
	while i < min(stop+1,len_cap):
		if i<start:
			_,_=cap.read()
		else:
			try:
				_,frame = cap.read()
				if r:
					frame = ndimage.rotate(frame, angle)
				if cr:
					frame=frame[crop_y0:crop_y1,crop_x0:crop_x1]
				#add scale?
				images.append(frame)
			except:
				print('Bad frame', i)
		i+=1
	cap.release()

	global frames
	frames = [np.array(x) for x in images]

	global fig, ax
	fig = plt.figure()
	ax = fig.gca()
	ax.imshow(frames[0], cmap= 'gray')
	ax.axis('off')
	slider_ax = plt.axes([0.15, 0.05, 0.7, 0.025])

	#add button "save"

	#For deprecated versions
	global slider
	try:
		slider = Slider(slider_ax, 'Frame', start, stop, valinit=0, valstep=1)
	except:
		slider = Slider(slider_ax, 'Frame', start, stop, valinit=0)

	slider.on_changed(update_slider)
	plt.show()

def update_slider(val):
	n = int(slider.val)
	ax.imshow(frames[n-st], cmap= 'gray')
	fig.canvas.draw_idle()
