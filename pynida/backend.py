#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy
from scipy import optimize, ndimage
import matplotlib.pyplot as plt
import matplotlib.patches

# For areas segmentation only
from skimage import measure  # <0.16
from skimage.segmentation import random_walker

import cv2

from pynida.simple_functions import thr_by_hist, basic_plot, thr_by_level, borderCheck, linear, find_minmax, \
    prepare_contours, \
    verify_contact, cone_volume
from pynida.classes.line import Line
from pynida.classes.point import Point


def oneshot_contours(img, savepath, plt_all=False):
    thr = thr_by_hist(img, 60, 10)  # img, default thr,sensitivity
    fin_contours = prepare_contours(img, thr)
    try_contours = verify_contact(fin_contours)
    results = {}
    if plt_all:
        plt.figure(frameon=False)
        plt.axis('off')
        plt.tight_layout()
        plt.subplots_adjust(right=0.98, left=0.02, top=0.98, bottom=0.02)
        plt.imshow(img, cmap='gray', interpolation='nearest')
        # plt.scatter(ml[0],ml[1],color= 'red')
        plt.savefig(savepath + '_c_initial.png')
        plt.close()

    if try_contours:
        aaa = fin_contours[0].T
        aaa = aaa[1][0]
        bbb = fin_contours[1].T
        bbb = bbb[1][0]
        if sum(aaa)/len(aaa) < sum(bbb)/len(bbb):
            top_c = fin_contours[0]
            bot_c = fin_contours[1]
        else:
            top_c = fin_contours[1]
            bot_c = fin_contours[0]
        '''
		if plt_all:
			plt.figure(frameon=False)
			plt.imshow(cr_test2, cmap= 'gray')
			cv2.drawContours(cr_test2,top_c, -1, (255,255,255), 3)
			cv2.drawContours(cr_test2,bot_c, -1, (255,255,255), 3)
			plt.axis('off')
			plt.tight_layout()
			plt.savefig(outdir+'/cont1/'+vhod[:-4]+'_'+make_str(frame_number)+'_m1.png')
			#plt.show()
			plt.close()
#		'''

        a, b, top_c = find_minmax(top_c, p='max')
        a2, b2, bot_c = find_minmax(bot_c, p='min')

        ts1 = min(a, b, key=(lambda p: p.x))
        sp1 = max(a, b, key=(lambda p: p.x))
        ts2 = min(a2, b2, key=(lambda p: p.x))
        sp2 = max(a2, b2, key=(lambda p: p.x))

        mp1 = Point.get_middle_point(ts1, ts2)
        mp2 = Point.get_middle_point(sp1, sp2)

        # TODO: rename it!
        line = Line(mp1, mp2)

        # TODO: rename it!
        def calculate_volumes(distances_from_contours_to_line, projected_contours):
            starting_point = projected_contours[0]
            dist = projected_contours.get_distances_to(starting_point)

            dl = np.diff(dist)
            radius = distances_from_contours_to_line[:-1]

            vol = np.sum(np.pi * dl * radius**2)
            invvol = np.sum(dl / radius**2 / np.pi)  # wout cones
            return vol, invvol


        top_c_l = line.orthogonal_projection(top_c)
        bot_c_l = line.orthogonal_projection(bot_c)

        dist_top = line.distances_to_points(top_c)
        dist_bot = line.distances_to_points(bot_c)

        vol_t, invvol_t = calculate_volumes(dist_top, top_c_l)
        vol_b, invvol_b = calculate_volumes(dist_bot, bot_c_l)

        # TODO: refactor volume corrections
        corr_vol_lb, corr_vol_lt = volume_corrections(bot_c_l[0], top_c_l[0], dist_bot[0],
                                                      dist_top[0], line, mp1, ts1, ts2, is_left=True)
        corr_vol_rb, corr_vol_rt = volume_corrections(bot_c_l[-1], top_c_l[-1], dist_bot[-1],
                                                      dist_top[-1], line, mp2, sp1, sp2, is_left=False)
        corrvol_t = corr_vol_lt + corr_vol_rt
        corrvol_b = corr_vol_lb + corr_vol_rb
        print("Vol ", vol_t, vol_b, corrvol_t, corrvol_b)

        btsc = Point.get_distance_between(ts1, ts2)
        bspc = Point.get_distance_between(sp1, sp2)
        mlc = Point.get_distance_between(mp1, mp2)

        dmx = abs(mp1.x - mp2.x)
        dmy = abs(mp1.y - mp2.y)
        # '''
        if plt_all:
            plt.figure(frameon=False)
            plt.axis('off')
            plt.tight_layout()
            plt.subplots_adjust(right=0.98, left=0.02, top=0.98, bottom=0.02)
            plt.imshow(img, cmap='gray')
            # plt.imshow(labels_fin, cmap= 'viridis', interpolation= 'nearest',alpha=.5)
            plt.scatter(top_c_l.x, top_c_l.y, color='red', marker='.')
            plt.scatter(top_c.x, top_c.y, color='blue', marker='.')
            plt.scatter(bot_c.x, bot_c.y, color='blue', marker='.')

            plt.scatter((ts1.x, ts2.x), (ts1.y, ts2.y), color='red', marker='x')
            plt.scatter((sp1.x, sp2.x), (sp1.y, sp2.y), color='red', marker='x')
            plt.scatter((mp1.x, mp2.x), (mp1.y, mp2.y), color='yellow', marker='v')
            # plt.scatter(ml[0],ml[1],color = 'red')
            plt.savefig(savepath + '_contours2.png')
            plt.close()
    # '''
    else:
        btsc = 0
        bspc = 0
        mlc = 0
        vol_t = 0
        invvol_t = 0
        corrvol_t = 0
        vol_b = 0
        invvol_b = 0
        corrvol_b = 0
        dmx = 0
        dmy = 0
    results['mid_x'] = dmx
    results['mid_y'] = dmy
    results['ts_len'] = btsc
    results['invvol_t'] = invvol_t
    results['sp_len'] = bspc
    results['mid_len'] = mlc
    results['vol_t'] = vol_t
    results['vol_b'] = vol_b
    results['invvol_b'] = invvol_b
    results['corrvol_b'] = corrvol_b
    results['corrvol_t'] = corrvol_t
    return results


# TODO: refactor volume corrections
def volume_corrections(bot_c_l: Point, top_c_l: Point, dist_bot: Point, dist_top: Point,
                       line, mp1: Point, ts1: Point, ts2: Point, is_left: bool):
    if ts1.x == ts2.x:
        # no correction needed
        return 0, 0

    if line.is_horizontal():
        ss = 1 if top_c_l.x > mp1.x else -1
        ss = ss if is_left else -ss
    else:
        ts_contact_line = Line(ts1, ts2)
        cosine = Line.cosine_distance(ts_contact_line, line)
        if cosine == 0:
            ss = 0
        elif cosine < 0:
            ss = -1
        else:
            ss = 1

    corr_vol_lt = ss*cone_volume(radius=dist_top,
                                 height=Point.get_distance_between(mp1, top_c_l))
    corr_vol_lb = -ss*cone_volume(radius=dist_bot,
                                  height=Point.get_distance_between(mp1, bot_c_l))
    return corr_vol_lb, corr_vol_lt


def oneshot_areas_proc(labels_fin, cr_test, output, plt_all=False):
    results = {}
    tip_label = 5

    # if labels_fin.find(1)>=0 and labels_fin.find(2)>=0:
    pod_label = 2
    sh_label = 1

    borders_p = ndimage.generic_filter(labels_fin, borderCheck, footprint=np.ones((3, 3)),
                                       extra_arguments=(-1, pod_label))
    i_borders_p = borders_p.nonzero()
    try:
        (a, _), _ = scipy.optimize.curve_fit(linear, i_borders_p[0], i_borders_p[1])
        p_angle = np.arctan(a)
    except:
        a = 0
        p_angle = 0

    results['substr_angle'] = p_angle

    borders_ts = ndimage.generic_filter(labels_fin, borderCheck, size=3, extra_arguments=(sh_label, tip_label))
    b_ts = np.nonzero(borders_ts)

    # Tip-sample border as a line from top to bottom points of the contour
    b_ts_c = np.array([[[x, y]] for x, y in zip(b_ts[0], b_ts[1])])
    if len(b_ts[0]) > 2:
        a_ts = np.array(b_ts_c[b_ts_c[:, :, 0].argmin()][0])
        b_ts = np.array(b_ts_c[b_ts_c[:, :, 0].argmax()][0])
        tip_sample_contact_length = np.sqrt(sum((b_ts - a_ts)**2))
    else:
        tip_sample_contact_length = 0
    # Tip-sample border as a line between most distant points of the contour (slow)
    '''
	if len(b_ts[0])>2:
		q = scipy.transpose(b_ts)
		a = scipy.spatial.distance.pdist(q)
		d = scipy.spatial.distance.squareform(a)
		tip_sample_contact_length = np.ndarray.max(d)
		ll=np.where(d == tip_sample_contact_length)
		print(ll)
		if ll:
			ll = ll[0]
		a = ll[0]
		b = ll[1]
		a_ts = q[a]
		b_ts = q[b]
		print(q[a],q[b])
	else:
		tip_sample_contact_length = 0
	'''
    print("tip_sample_contact_length ", tip_sample_contact_length)
    results["tip_sample_contact_length"] = tip_sample_contact_length

    # edge_ts = ndimage.generic_filter(labels_fin,edgeCheck,size = 3,extra_arguments = (sh_label,tip_label))
    # i_ets = edge_ts.nonzero()
    # print(i_ets)

    borders_ps = ndimage.generic_filter(labels_fin, borderCheck, footprint=np.ones((3, 3)),
                                        extra_arguments=(sh_label, pod_label))
    print("Borders ps!")
    b_ps = np.nonzero(borders_ps)

    # Surface-sample border as a line from top to bottom points of the contour
    b_ps_c = np.array([[[x, y]] for x, y in zip(b_ps[0], b_ps[1])])
    if len(b_ps[0]) > 2:
        a_ps = np.array(b_ps_c[b_ps_c[:, :, 0].argmin()][0])
        b_ps = np.array(b_ps_c[b_ps_c[:, :, 0].argmax()][0])
        sample_surface_contact_length = np.sqrt(sum((b_ps - a_ps)**2))
    else:
        sample_surface_contact_length = 0
    # Surface-sample border as a line between most distant points of the contour (slow)
    # print(b_ts)
    '''
	if len(b_ps[0])>2:
		q = scipy.transpose(b_ps)
		a = scipy.spatial.distance.pdist(q)
		d = scipy.spatial.distance.squareform(a)
		sample_surface_contact_length = np.ndarray.max(d)
		ll=np.where(d == sample_surface_contact_length)
		print(ll)
		if ll:
			ll = ll[0]
		a = ll[0]
		b = ll[1]
		a_ps = q[a]
		b_ps = q[b]
		print(q[a],q[b])
	else:
		sample_surface_contact_length = 0
	'''
    print("sample_surface_contact_length ", sample_surface_contact_length)
    results["sample_surface_contact_length"] = sample_surface_contact_length

    len_mid = 0
    if (sample_surface_contact_length > 0.) and (tip_sample_contact_length > 0.):
        print("Both contacts are present")
        mpl = (a_ts + b_ts)/2.
        mpr = (a_ps + b_ps)/2.
        mp = mpr - mpl
        dx_tip_sample = abs(a_ts[1] - b_ts[1])
        dx_sample_surface = abs(a_ps[1] - b_ps[1])
        len_mid = np.sqrt(sum(mp**2))
        print("len_mid", len_mid)
        results["mid_projection_x"] = abs(mp[1])
        results["mid_projection_y"] = abs(mp[0])
        results["len_mid"] = len_mid
        results["dx_tip_sample"] = dx_tip_sample
        results["dx_sample_surface"] = dx_sample_surface
        if plt_all:
            _, ax = plt.subplots()
            ax.imshow(cr_test, cmap=plt.cm.gray, alpha=0.5)
            ax.imshow(labels_fin, cmap='viridis', interpolation='nearest', alpha=.5)
            plt.subplots_adjust(right=0.98, left=0.02, top=0.98, bottom=0.02)
            plt.plot((a_ts[1], b_ts[1]), (a_ts[0], b_ts[0]), color='red', marker='.')
            plt.plot((a_ps[1], b_ps[1]), (a_ps[0], b_ps[0]), color='red', marker='.')
            plt.axis('off')
            plt.savefig(output + '_Test_area_fin.png')
            plt.close()
    #	D_mlen.append(len_mid)
    else:
        results["mid_projection_x"] = 0
        results["mid_projection_y"] = 0
        results["len_mid"] = 0
        results["dx_tip_sample"] = 0
        results["dx_sample_surface"] = 0

    sph = cr_test.copy()
    sph[labels_fin != sh_label] = 0
    sph2 = sph.copy()
    sph[labels_fin == sh_label] = 1

    #######
    # http://scikit-image.org/docs/dev/auto_examples/segmentation/plot_regionprops.html#sphx-glr-auto-examples-segmentation-plot-regionprops-py

    lll = measure.label(sph)
    # plt.imshow(lll)
    # plt.show()
    # plt.close()
    try:
        prop = measure.regionprops(lll, coordinates='xy')
    except:
        prop = measure.regionprops(lll)
    if prop:
        prop = prop[0]
        y0, x0 = prop.centroid
        orientation = prop.orientation
        minr, minc, maxr, maxc = prop.bbox
        bx = (minc, maxc, maxc, minc, minc)
        by = (minr, minr, maxr, maxr, minr)

        square_x = maxr - minr
        square_y = maxc - minc
        ellipse_major = prop.major_axis_length
        ellipse_minor = prop.minor_axis_length
        area = prop.area
        if plt_all:
            _, ax = plt.subplots()
            ax.imshow(sph2, cmap=plt.cm.gray)

            loc_x1 = x0 + np.cos(orientation)*0.5*prop.major_axis_length
            loc_y1 = y0 - np.sin(orientation)*0.5*prop.major_axis_length
            loc_x2 = x0 - np.sin(orientation)*0.5*prop.minor_axis_length
            loc_y2 = y0 - np.cos(orientation)*0.5*prop.minor_axis_length

            ax.plot((x0, loc_x1), (y0, loc_y1), '-r', linewidth=1.5)
            ax.plot((x0, loc_x2), (y0, loc_y2), '-r', linewidth=1.5)
            ax.plot(x0, y0, '.g', markersize=15)
            ax.add_artist(matplotlib.patches.Ellipse(xy=(x0, y0), width=ellipse_major, height=ellipse_minor,
                                                     angle=-orientation*180./np.pi, alpha=.5, fc=None, ec='red',
                                                     lw=2))

            ax.plot(bx, by, '-b', linewidth=1)
            # ax.axis((0,2048,0,1600))
            # plt.show()
            plt.axis('off')
            # plt.tight_layout()
            plt.subplots_adjust(right=0.98, left=0.02, top=0.98, bottom=0.02)
            plt.savefig(output + '_Test_area_particle.png')
            plt.close()
    # '''
    # print(orientation)
    else:
        square_x = 0
        square_y = 0
        ellipse_major = 0
        ellipse_minor = 0
        orientation = 0
        area = 0

    results["square_x"] = square_x
    results["square_y"] = square_y
    results["ellipse_major"] = ellipse_major
    results["ellipse_minor"] = ellipse_minor
    results["orientation"] = orientation
    results["area"] = area

    return results


def oneshot_areas_select(img, output_dir, fr_number,
                         borders_remove1, borders_remove2,
                         right_cut2=0., plt_all=False,
                         blur_val=3, sm_blur_val=21,
                         min_area1=1000, start_coeff1=0.5, step_coeff1=0.01,
                         min_area2=500, start_coeff2=0.5, step_coeff2=0.001):
    cr_test = cv2.GaussianBlur(img, (blur_val, blur_val), 0)
    #	smoothest = cv2.GaussianBlur(cr_test,(sm_blur_val,sm_blur_val),0)

    first_thr = thr_by_hist(cr_test, 60, 10)  # img, default thr,sensitivity

    #	tmp_thr = 100
    #	tmp_thr2 = 100
    #	second_thr=100#For markers
    #	tmp_thr2=tmp_thr
    #	second_thr = tmp_thr

    #	if override_first_thr:
    #		first_thr=150
    #		tmp_thr = 150

    shpod = cr_test.copy()
    shpod[shpod < first_thr] = 0
    shpod2 = shpod.copy()

    if plt_all:
        basic_plot(cr_test, f'{output_dir}_{fr_number}_area_initial.png')

    shpod_tmp = shpod.copy()
    shpod_tmp = cv2.GaussianBlur(shpod, (blur_val, blur_val), 0)
    # shpod_tmp[shpod_tmp<tmp_thr]=0
    '''
	if borders_remove:
	order: 'l','r','b','t'
	'''
    if borders_remove1[0]:  # _l:
        shpod_tmp[:, 0] = 0
    if borders_remove1[1]:  # _r:
        shpod_tmp[:, -1] = 0
    if borders_remove1[2]:  # _b:
        shpod_tmp[-1, :] = 0
    if borders_remove1[3]:  # _t:
        shpod_tmp[0, :] = 0

    sh_distance = ndimage.distance_transform_edt(shpod_tmp)

    if plt_all:
        basic_plot(sh_distance, f'{output_dir}_{fr_number}_area_dist1.png', cmap='viridis')

    # size of area!
    thr_dist = thr_by_level(sh_distance, min_area1,
                            step_coeff1, start_coeff1)  # Threshold level for dist map. min area in px, step, init_thr

    shpod_tmp[sh_distance < thr_dist] = 0
    shpod_tmp[sh_distance >= thr_dist] = 1

    markers = measure.label(shpod_tmp.T)
    markers = markers.T
    markers[markers >= 3] = 0
    markers[cr_test < first_thr] = -1

    if plt_all:
        basic_plot(markers, f'{output_dir}_{fr_number}_area_marker1.png', cmap='viridis')

    try:
        labels_shpod = random_walker(shpod, markers, beta=10, mode='bf')
    except:
        labels_shpod = markers
        print("Error in 1st random walker")

    if plt_all:
        basic_plot(labels_shpod, f'{output_dir}_{fr_number}_area_rw1.png', cmap='viridis')

    #	tmp_thr2=100#for segmentation

    shpod2[labels_shpod == 1] = 0
    shpod_tmp2 = shpod2.copy()
    shpod_tmp2 = cv2.blur(shpod_tmp2, (sm_blur_val, sm_blur_val))
    '''
	order: 'l','r','b','t'
	'''
    if borders_remove2[0]:  # _l:
        shpod_tmp2[:, 0] = 0
    if borders_remove2[1]:  # _r:
        shpod_tmp2[:, -1] = 0
    if borders_remove2[2]:  # _b:
        shpod_tmp2[-1, :] = 0
    if borders_remove2[3]:  # _t:
        shpod_tmp2[0, :] = 0

    if right_cut2 > 0.:
        shpod_tmp2[:, -int(right_cut2/100*len(shpod_tmp2[0])):] = 0

    #	shpod_tmp2[cr_test < second_thr]=0

    if plt_all:
        basic_plot(shpod2, f'{output_dir}_{fr_number}_area_init2.png')

    #	shpod_tmp2[shpod2<tmp_thr2]=0
    #	shpod_tmp2[shpod2>=tmp_thr2]=1
    #	shpod_tmp2 = cv2.blur(shpod_tmp2,(blur_val,blur_val))

    sh_distance2 = ndimage.distance_transform_edt(shpod_tmp2)

    if plt_all:
        basic_plot(sh_distance2, f'{output_dir}_{fr_number}_area_dist2.png', cmap='viridis')

    tmp_thr2 = thr_by_level(sh_distance2, min_area2, step_coeff2,
                            start_coeff2)  # Threshold level for dist map. min area in px

    shpod_tmp2 = (sh_distance2 >= tmp_thr2)

    markers2 = measure.label(shpod_tmp2.T)

    markers2 = markers2.T
    markers2[markers2 >= 3] = 0
    markers2[cr_test < first_thr] = -1
    markers2[labels_shpod == 1] = -1

    if plt_all:
        basic_plot(markers2, f'{output_dir}_{fr_number}_area_markers2.png', cmap='viridis')

    try:
        labels_shpod2 = random_walker(shpod2, markers2, beta=100, mode='bf')
    except:
        labels_shpod2 = markers2
        print("Error in 2nd random walker")
    # labels_fin = labels_shpod2.copy()
    labels_shpod2[labels_shpod == 1] = 5
    # tip_label = 5

    if plt_all:
        basic_plot(labels_shpod2, f'{output_dir}_{fr_number}_area_rw2.png', cmap='viridis')
    # plt.scatter(ml[0],ml[1],color= 'red')

    return labels_shpod2
