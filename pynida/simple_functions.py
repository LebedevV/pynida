#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy
import cv2
import matplotlib.pyplot as plt
from skimage import measure #<0.16
from typing import Tuple

import pynida.classes
from pynida.classes.point import Point
from pynida.classes.points import Points


def get_n_frame(f, start, stop, n):

	n = max(0, n, start)
	n = min(n, stop)

	cap = cv2.VideoCapture(f)

	for i in range(stop):
		img = cap.read()[1]
		if i >= n:
			break

	cap.release()
	return img

def make_str(n):
	n = str(n)
	while len(n) < 5:
		n = '0' + n
	return n

def sneddon(x, x0, y0, alpha):
	x = x - x0
	x[x < 0] = 0
	y = y0 + alpha * np.power(x, 1.5)
	return y

def linear(x, a, b): #for linear fit
	y = a * x + b
	return y

def thr_by_hist(img, threshold, edge):
	hist = cv2.calcHist([img], [0], None, [256], [0, 256])
	hist = np.squeeze(hist)[:-1]  # reject max value at last bin

	hist_diff = -np.diff(hist)
	peak_position = np.argmax(hist)
	# find the rightmost point of the background peak
	for i in range(peak_position, len(hist) - 2):
		if hist_diff[i] > edge > hist_diff[i+1]:
			return i + 10
	print('Something wrong with threshold')  # error msg
	print(f'Basic threshold {threshold}')
	return threshold


def thr_by_level(img, minimal_area, step, constant):
	thr = cv2.minMaxLoc(img)[1] #min_val,max_val,min_loc,max_loc
	while 1:
		constant -= step
		tthr = thr * constant
		if constant <= 0.01:
			print('Error - low area threshold level')
			break

		mask = (img >= tthr)
		m = measure.label(mask, background=0)
		if min(np.count_nonzero(m == 1), np.count_nonzero(m == 2)) > minimal_area:
			print("Threshold done well")
			break
	return tthr


def basic_plot(fr, n, cmap='gray'):
	plt.figure(frameon=False)
	plt.axis('off')
	plt.tight_layout()
	plt.subplots_adjust(right=0.98, left=0.02, top=0.98, bottom=0.02)
	plt.imshow(fr, cmap=cmap, interpolation='nearest')
	#plt.scatter(ml[0],ml[1],color= 'red')
	plt.savefig(n)
	plt.close()


def find_minmax(c, p='max', xcutoff=100) -> Tuple[Point, Point, Points]: #Add minmax control, unite with find_max.

	ini = c.copy()
	#First extremum
	if p == 'min':
		most1 = Point.from_tuple(tuple(c[c[:, :, 1].argmin()][0]))
		i_res = c[:,:,1].argmax()
		thr = 10000
	elif p == 'max':
		most1 = Point.from_tuple(tuple(c[c[:, :, 1].argmax()][0]))
		i_res = c[:,:,1].argmin()
		thr = 0
	#print("most1",most1)

	#Distance to the first extremum
	dist = abs(c[:,:,0] - most1.x)[:, 0]

	#removal of the first extremum
	c[:,:,1][dist < xcutoff] = thr #cut out the area around firs extremum

	#second extremum
	if p == 'min':
		most2 = Point.from_tuple(tuple(c[c[:, :, 1].argmin()][0]))
	elif p == 'max':
		most2 = Point.from_tuple(tuple(c[c[:, :, 1].argmax()][0]))

	#print("most2",most2)

	ini2 = ini.copy()
	ini2 = ini2.tolist()
	ini2 = ini2[i_res:] + ini2[:i_res] #outer edge removal

	#point-by-point comparison and adding to the final contour
	if (most1.x < most2.x and p == 'min') or (most1.x > most2.x and p == 'max'):
		ini2.reverse()
	add = False
	fin = Points()

	for t in ini2:
		p = Point.from_tuple(t[0])
		if p == most1:
			print("P1_top found!")
			add = True
		if add:
			fin.append(Point(p.x, p.y))
		if p == most2:
			print("P2_top found!")
			add = False

	return most1, most2, fin



def find_minmax_test(c, p='max', xcutoff=50): #Add minmax control, unite with find_max.
	ini = c.copy()
	#First extremum
	if p == 'min':
		peak_x = c[:,:,0].argmin()
		#print('peak_x',peak_x)
		most1 = tuple(c[peak_x][0])
		i_res = c[:,:,0].argmax()
		#print('i',i_res)
		thr = 10000
	elif p == 'max':
		peak_x = c[:,:,0].argmax()
		#print('peak_x',peak_x)
		most1 = tuple(c[peak_x][0])
		i_res = c[:,:,0].argmin()
		thr = 0
		#print('i',i_res)
	#print("most1",most1)

#	'''
	#Distance to the first extremum
	dist = abs(c[:,:,0] - most1[0])[:, 0]
	#removal of the first extremum
	c[:,:,1][dist < xcutoff] = thr #cut out the area around firs extremum
#	'''

	'''
	c = c.tolist()
	c = c[:peak_x-xcutoff] + c[peak_x+xcutoff:]
	#c = np.concatenate(c[:,:,:peak_x-xcutoff],c[:,:,peak_x+xcutoff:])
	c = np.array(c)
	print('c',c)
#	'''

	#second extremum
	if p == 'min':
		most2 = tuple(c[c[:,:,0].argmin()][0])
	elif p == 'max':
		most2 = tuple(c[c[:,:,0].argmax()][0])

	#print("most2",most2)

#	'''
	cut = ini.copy()
	cut = cut.tolist()
	#cut = cut[i_res:] + cut[:i_res] #outer edge removal
	#point-by-point comparison and adding to the final contour
#	if (most1[0]<most2[0] and p == 'min') or (most1[0]>most2[0] and p == 'max'):
#		cut.reverse()
	a = cut.index([list(most1)])
	b = cut.index([list(most2)])
	if p == 'min':
		cut = cut[a:b]
	else:
		cut = cut[a:] + cut[:b]
#	cut = cut[min(a,b):max(a,b)]
	print('cut after', cut)
#	'''
	cut = np.array(cut)
	fin_x = cut[:,:,0]
	fin_y = cut[:,:,1]

	return most1[0], most1[1], most2[0], most2[1], np.array(fin_x), np.array(fin_y)

def parse_video_data(f, sep, bckgr):
	with open(f) as file:
		data = file.read()
	data = data.split('\n')
	info = data[:data.index(sep)]
	data = data[data.index(sep)+2:-1]
	time = []
	disp = []
	for row in data:
		row = row[:-2]
		temp = row.split('\t')
		#print(temp)
		temp = np.array(temp, dtype= 'float')
		time.append(temp[0]) # frame numbers
		if bckgr:
			disp.append(temp[1] - temp[4]) # displacement in X axis
		else:
			disp.append(temp[1])
	pnmi= '###'
	for s in info:
		if s.startswith('pnm'):
			pnmi = s.split('\t')[1]
	return time, disp, pnmi


def parse_old_data(filename):
	with open(filename, encoding='ISO-8859-1') as file:
		data = file.readlines()
		
	t = 0
	raw = False
	for d in data:
		if 'Depth (nm)	Load (µN)' in d:
			t = data.index(d) + 1
			raw = True
		if '---' in d:
			t = data.index(d) + 2
	data = data[t:]
	time = []
	x = []
	f = []
	for row in data:
		ss = row.split('\t')
		ss = np.array(ss, dtype='float')
		if raw:
			x.append(ss[0]) # depth, nm
			f.append(ss[1]) # load, microN
			time.append(ss[2]) # time, s
		else:
			x.append(ss[1]) # depth, nm
			f.append(ss[2]) # load, microN
			time.append(ss[0]) # time, s
	return time, x, f


def parse_hys_data(filename):
	with open(filename, encoding='ISO-8859-1') as file:
		data = file.readlines()
	data = data[7:]
	time = []
	x = []
	f = []
	for row in data:
		ss = row.split('\t')
		ss = np.array(ss, dtype='float')
		x.append(ss[0]) # depth, nm
		f.append(ss[1]) # load, microN
		time.append(ss[2]) # time, s
	return time, x, f

def parse_corr_data(filename, separator):
	with open(filename, encoding='ISO-8859-1') as file:
		data = file.readlines()
	data = data[data.index(separator + '\n') + 2: -1]
	time = []
	x = []
	f = []
	for row in data:
		ss = row.split('\t')
		ss = np.array(ss, dtype='float')
		x.append(ss[1]) # depth, nm
		f.append(ss[2]) # load, microN
		time.append(ss[0]) # time, s
	return time, x, f

def interpolate_local(x, y, N):
	if len(x) != len(y):
		raise NameError('X and Y have different length')
	f = scipy.interpolate.interp1d(x, y, 'cubic')
	x_interpolated = np.linspace(min(x), max(x), N)
	y_interpolated = f(x_interpolated)
	return x_interpolated, y_interpolated

def Squeeze(array):
	span = float(max(array) - min(array))
	array_new = [x / span  for x in array]
	return array_new

def Convert(hys_time, hys_disp, vis_time, vis_disp, k): # 1 - hys, 2 - vis
	hys_time = np.array(hys_time)
	hys_disp = np.array(hys_disp)
	vis_time = np.array(vis_time)
	vis_disp = np.array(vis_disp)

	hys_disp_new = (k['disp_shift'] + hys_disp * k['disp_scale'] + hys_time * k['line_scale'] + 
			 k['sine_amplitude'] * np.sin(hys_time / k['sine_frequency'] + k['sine_shift']))
	vis_time_new = (vis_time + k['time_shift']) * k['time_scale']
	return hys_time, hys_disp_new, vis_time_new, vis_disp

def ErrorFunc(k, hys_time, hys_disp, vis_time, vis_disp, func):
	N = len(hys_time)
	if N != len(vis_time) or N != len(hys_disp) or N != len(vis_disp):
		raise NameError('X and Y have different length')
	hys_time, hys_disp_new, vis_time_new, vis_disp = func(hys_time, hys_disp, vis_time, vis_disp, k)
	for array in [hys_time, hys_disp_new, vis_time_new, vis_disp]:
		array = Squeeze(array)
	custom_error = (hys_time - vis_time_new) ** 2 + (hys_disp_new - vis_disp) ** 2
	custom_error = custom_error.sum() / N
	return custom_error

def eab(ea, eb):
	x = 1./ea + 1./eb
	return 1./x

def prepare_contours(init_frame, thr):
	_, tr = cv2.threshold(init_frame, thr, 255, cv2.THRESH_BINARY_INV)
	contours, _ = cv2.findContours(tr, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
	print("N contours= ", len(contours))

	#Contours sorted by len
	lens = [len(i) for i in contours]
	fin_c = [y for _, y in sorted(list(zip(lens, contours)), key = lambda x: x[0], reverse=True)]

	return fin_c

def verify_contact(fcs):
	r = False
	print("Verify contacts")
	if len(fcs) > 1:
		print(len(fcs[0]))
		print(len(fcs[1]))
		#print(len(fcs[2][0]))
		if len(fcs[0]) / len(fcs[1]) <= 10.:
			 r = True
	return r


def borderCheck(x,a,b):
	#c= (len(x)-1)/2		#Central point of x array
	c = 4
	if (a in x) and (x[c] == b):	#are there any pixels with "a"
#		print("Border!")
		y = 1#np.sum(x)/np.count_nonzero(x)
	else:
		y = 0
	return y


def cone_volume(radius, height):
    area = np.pi*radius**2
    return area*height/3
