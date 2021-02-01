import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy import optimize

from pynida.simple_functions import interpolate_local, Convert, ErrorFunc

# +1. create objects/dicts hys data and video data with attributes time and displacement
# +2. move sliders into sliders attribute (dict)
# (later)3. import **args - not dictionary
# +4. solve duplicate slider reading in update and reset
# +5. change parameters from array to dictionary (involves changing Convert)
# 6. change border_idx and borders from array to dictionary

class DynamicObject():
    def keys(self):
        return [x for x in dir(self) if not x.startswith('__') and x != "keys"]
    pass

class WindowComparison():
    def __init__(self, arguments):
        self.data = DynamicObject()
        self.data.video = DynamicObject()
        self.data.hys = DynamicObject()
        self.sliders = DynamicObject()
        self.sliders.parameters = DynamicObject()
        self.sliders.points = DynamicObject()
        self.plots = DynamicObject()

        self.data.video.time = arguments['time_v']
        self.data.video.new_time = []
        self.data.video.disp = arguments['disp_v']
        self.data.hys.time = arguments['time_h']
        self.data.hys.disp = arguments['disp_h']
        self.data.hys.new_disp = []

        self.pnmi = arguments['pnmi']

        self.vlb = None
        self.vrb = None
        self.hlb = None
        self.hrb = None
        self.border_idx = None
        self.parameters = {}
        self.parameters_old = {}
        self.slider_parameter_pairs = ()

    def get_parameters(self):
        return self.data.hys.time, self.data.hys.new_disp, self.data.video.new_time, self.data.video.disp, self.parameters

    def run(self):
        # draw window

        # Create plt window
        _, ax = plt.subplots()
        plt.xlabel("$Time, s$")
        plt.ylabel("$Displacement, nm$")
        indent_plot = 0.55
        plt.subplots_adjust(left=0.1, bottom = indent_plot,top=.98)

        # self.data.video.time is a frame number; it has to be the same order of value as self.data.hys.time
        scale = max(self.data.hys.time) - min(self.data.hys.time)
        scale /= (max(self.data.video.time) - min(self.data.video.time))
        #self.data.video.time = [x*scale for x in self.data.video.time]
        # plotting of 2 curves
        self.plots.video, = plt.plot(np.array(self.data.video.time) * scale, self.data.video.disp)
        self.plots.hys, = plt.plot(self.data.hys.time, self.data.hys.disp)
        # plot endpoints
        # fraction of points excluded from each side; slider values
        # borders = [0, 0, 0, 0]
        # indices of border points
        self.border_idx = [0, len(self.data.hys.time) - 1,
                                                   0, len(self.data.video.time) - 1]

        # plotting border points
        self.hlb, self.hrb = plt.plot(
            self.data.hys.time[self.border_idx[0]],
            self.data.hys.disp[self.border_idx[0]],
            'o',
            self.data.hys.time[self.border_idx[1]],
            self.data.hys.disp[self.border_idx[1]],
            'o',
            color = 'green')
        
        self.vlb, self.vrb = plt.plot(
            self.data.video.time[self.border_idx[2]] * scale,
            self.data.video.disp[self.border_idx[2]],
            'o',
            self.data.video.time[self.border_idx[3]] * scale,
            self.data.video.disp[self.border_idx[3]],
            'o',
            color = 'red')

        # axes for sliders
        axcolor = 'lightgoldenrodyellow'
        indent1 = indent_plot - 0.15 # indent1
        indent2 = 0.025 # indent2
        axshift  = plt.axes([0.1, indent1,                0.8, 0.01], facecolor=axcolor, title="Linear transform")
        axshifty = plt.axes([0.1, indent1 - indent2,      0.8, 0.01], facecolor=axcolor)
        axa      = plt.axes([0.1, indent1 - indent2 * 2,  0.8, 0.01], facecolor=axcolor)
        axscalex = plt.axes([0.1, indent1 - indent2 * 4,  0.8, 0.01], facecolor=axcolor, title="Scale factors")
        axscaley = plt.axes([0.1, indent1 - indent2 * 5,  0.8, 0.01], facecolor=axcolor)
        axsina   = plt.axes([0.1, indent1 - indent2 * 7,  0.8, 0.01], facecolor=axcolor, title="Sin(x) parameters")
        axsinf   = plt.axes([0.1, indent1 - indent2 * 8,  0.8, 0.01], facecolor=axcolor)
        axsinsh  = plt.axes([0.1, indent1 - indent2 * 9,  0.8, 0.01], facecolor=axcolor)
        axhyslb  = plt.axes([0.1, indent1 - indent2 * 11, 0.8, 0.01], facecolor=axcolor, title="Borders")
        axhysrb  = plt.axes([0.1, indent1 - indent2 * 12, 0.8, 0.01], facecolor=axcolor)
        axvislb  = plt.axes([0.1, indent1 - indent2 * 13, 0.8, 0.01], facecolor=axcolor)
        axvisrb  = plt.axes([0.1, indent1 - indent2 * 14, 0.8, 0.01], facecolor=axcolor)

        # Sliders
        x = self.sliders.parameters
        # shift and scale
        x.shiftx = Slider(axshift,  'ShiftX', -300., 300.0, valinit=0.)
        x.shifty = Slider(axshifty, 'ShiftY', -100., 100, valinit=0.)
        x.scalex = Slider(axscalex, 'ScaleX', 0.75 * scale, 5 * scale, valinit=scale)
        x.scaley = Slider(axscaley, 'ScaleY', 0.5, 5., valinit=.89)
        # drift function parameters
        x.a     = Slider(axa,         'k',   -5,   5, valinit=0.)
        x.sina  = Slider(axsina,    'Amp',    0, 500, valinit=0.)
        x.sinf  = Slider(axsinf,   'Freq',    1, 500, valinit=150.)
        x.sinsh = Slider(axsinsh, 'Phase', -50., 50., valinit=0.)
        # pair sliders with drift coefficients for future use
        self.slider_parameter_pairs = (
            (x.scalex, 'time_scale'),
            (x.shiftx, 'time_shift'),
            (x.scaley, 'disp_scale'),
            (x.shifty, 'disp_shift'),
            (x.a, 'line_scale'),
            (x.sina, 'sine_amplitude'),
            (x.sinf, 'sine_frequency'),
            (x.sinsh, 'sine_shift')
        )

        # border parameters
        x = self.sliders.points 
        x.hyslb = Slider(axhyslb,  'Left hys',   0, 50, valinit=0)
        x.hysrb = Slider(axhysrb, 'Right hys', -50,  0, valinit=0)
        x.vislb = Slider(axvislb,  'Left vis',   0, 50, valinit=0)
        x.visrb = Slider(axvisrb, 'Right vis', -50,  0, valinit=0)

        # Buttons
        axrefine = plt.axes([0.91,  0.9, 0.07, 0.06])
        axreset  = plt.axes([0.91, 0.83, 0.07, 0.06])
        self.brefine = Button(axrefine, 'refine')
        self.breset  = Button(axreset, 'reset')
        self.brefine.on_clicked(self.refine)
        self.breset.on_clicked(self.reset)

        # set slider fonts
        for ax in [ax,axshift,axshifty,axa,axscalex,axscaley,axsina,axsinf,axsinsh,axhyslb,axhysrb,axvislb,axvisrb]:
                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(10)

        # connect with update function
        for x in [self.sliders.parameters, self.sliders.points]:
                for slider in x.keys():
                        getattr(x, slider).on_changed(self.update)

        plt.show()

    def refine(self, event, N = 500):
        self.__update_parameters()
        self.__update_borders()

        loc_x1 = self.data.hys.time[self.border_idx[0]:self.border_idx[1]]
        loc_y1 = self.data.hys.disp[self.border_idx[0]:self.border_idx[1]]
        loc_x2 = self.data.video.time[self.border_idx[2]:self.border_idx[3]]
        loc_y2 = self.data.video.disp[self.border_idx[2]:self.border_idx[3]]

        x1_interpolated, y1_interpolated = interpolate_local(loc_x1, loc_y1, N)
        x2_interpolated, y2_interpolated = interpolate_local(loc_x2, loc_y2, N)

        solution = optimize.minimize(
                ErrorFunc, 
                list(self.parameters.values()),
                args = (x1_interpolated, 
                        y1_interpolated, 
                        x2_interpolated, 
                        y2_interpolated, 
                        self.__convert_wrapper),
                bounds = ((0.5, 20), (-150, 150), (0.1, 10), (-100, 100),
                            (-5, 5),    (0, 500),  (1, 500),   (-0, 50)))
        print(solution)
        self.parameters_old = self.parameters.copy()
        # convert result of optimization from list into dictionary
        for key, value in zip(self.parameters.keys(), solution.x):
                self.parameters[key] = value

        self.__update_data()
        self.__set_slider_values()
        self.__update_plot()
        plt.draw()

    def __convert_wrapper(self, hys_time, hys_disp, vis_time, vis_disp, k):
        '''Wrapper for Convert function: 
            converts array k into a dictionary and calls Convert
        '''
        parameters = {}
        for i, pair in enumerate(self.slider_parameter_pairs):
            key = pair[1]
            parameters[key] = k[i]
        return Convert(hys_time, hys_disp, vis_time, vis_disp, parameters)

    def reset(self, event):
        self.parameters = self.parameters_old
        self.__set_slider_values()
        self.__update_data()
        self.__update_plot()
        plt.draw()

    def update(self, val):
        self.__update_parameters()
        self.__update_borders()
        self.__update_data()
        self.__update_plot()
        plt.draw()


    def __update_data(self):
        _, self.data.hys.new_disp, self.data.video.new_time, _ = \
                Convert(self.data.hys.time, 
                        self.data.hys.disp, 
                        self.data.video.time, 
                        self.data.video.disp, 
                        self.parameters)

    def __update_borders(self):
        x = self.sliders.points
        values = [x.hyslb.val, x.hysrb.val, x.vislb.val, x.visrb.val]
        borders = [float(x / 100) for x in values]
        self.border_idx = [
            len(self.data.hys.time) * borders[0],
            len(self.data.hys.time) * (1 + borders[1]) - 1,
            len(self.data.video.time) * borders[2],
            len(self.data.video.time) * (1 + borders[3]) - 1
        ]
        self.border_idx = [int(x) for x in self.border_idx]

    def __update_parameters(self):
        for slider, key in self.slider_parameter_pairs:
            self.parameters[key] = float(slider.val)

    def __set_slider_values(self):
        for slider, key in self.slider_parameter_pairs:
            slider.set_val(self.parameters[key])

    def __update_plot(self):
        self.vlb.set_xdata(self.data.video.new_time[self.border_idx[2]])
        self.vrb.set_xdata(self.data.video.new_time[self.border_idx[3]])
        self.vlb.set_ydata(self.data.video.disp[self.border_idx[2]])
        self.vrb.set_ydata(self.data.video.disp[self.border_idx[3]])
        self.hlb.set_ydata(self.data.hys.new_disp[self.border_idx[0]])
        self.hrb.set_ydata(self.data.hys.new_disp[self.border_idx[1]])
        self.hlb.set_xdata(self.data.hys.time[self.border_idx[0]])
        self.hrb.set_xdata(self.data.hys.time[self.border_idx[1]])
        self.plots.video.set_xdata(self.data.video.new_time)
        self.plots.hys.set_ydata(self.data.hys.new_disp)
