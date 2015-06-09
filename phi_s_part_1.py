import subprocess
import math
import numpy
import sys
import os
import glob
import json
import ROOT
from xboa.algorithms.peak_finder import WindowPeakFinder
from xboa.algorithms.peak_finder import RefinePeakFinder
from xboa.algorithms.peak_finder import UphillDownhillPeakFinder
from xboa.algorithms.smoothing import GaussianSmoothing
import xboa.common as common

#TODO: get the peak
#TODO:

DT = 0.2

def load_data_from_listing(file_listing):
    data_list = []
    for filename in [file_listing["rf"], file_listing["ac"], file_listing["dc"]]:
        print filename
        fin = open(os.path.expandvars(filename))
        time_list, voltage_list = [], []
        for line in fin.readlines():
            words = line.split(',')
            words = [x for x in words if x != ""]
            try:
                time_list.append(float(words[0]))
                voltage_list.append(float(words[1]))
            except ValueError:
                 pass
        print "  read", len(time_list), "values"
        data = {}
        data["filename"] = os.path.split(filename)[-1]
        data["time_list"] = time_list
        data["voltage_list"] = voltage_list
        data["signal"] = file_listing["signal"]
        data["rf_voltage"] = file_listing["v_in"]
        data_list.append(data)
        print "Loaded", data["filename"], data["rf_voltage"], data["signal"]
    pair = {"RF":data_list[0], "signal":data_list[1], "dc_signal":data_list[2]}
    if pair["RF"]["rf_voltage"] != pair["signal"]["rf_voltage"] or \
       pair["RF"]["rf_voltage"] != pair["dc_signal"]["rf_voltage"]:
        raise ValueError("Failed to split data_list properly")
    return pair

def load_data(glob_name):
    data_list = []
    for filename in glob.glob(glob_name):
        fin = open(filename)
        time_list, voltage_list = [], []
        for line in fin.readlines():
            words = line.split(',')
            try:
                time_list.append(float(words[0]))
                voltage_list.append(float(words[1]))
            except ValueError:
                 pass
        data = {}
        data["filename"] = os.path.split(filename)[-1]
        data["time_list"] = time_list
        data["voltage_list"] = voltage_list
        data["signal"] = filename[filename.find("_63_")+4:filename.find("_ch")]
        data["rf_voltage"] = float(filename[filename.find("RF_")+3:filename.find("_foil")])
        data_list.append(data)
        print "Loaded", data["filename"], data["rf_voltage"], data["signal"]
    data_list = sorted(data_list, key = lambda data: [data["rf_voltage"], data["signal"]])
    data_list_pairs = []
    for i in range(0, len(data_list), 2):
        pair = {"RF":data_list[i], "signal":data_list[i+1]}
        if pair["RF"]["rf_voltage"] != pair["signal"]["rf_voltage"]:
            raise ValueError("Failed to split data_list properly")
        data_list_pairs.append(pair)
    return data_list_pairs

def mean_sigma_cut(data):
    sigma_list = data
    while True:
        n_sigma_old = len(sigma_list)
        mean = numpy.mean(sigma_list)
        sigma = numpy.std(sigma_list)
        sigma_list = [x for x in sigma_list if abs(mean-x) < 3.*sigma]
        if len(sigma_list) == n_sigma_old:
            break
    print "    Delta",
    print "data len:", len(data), len(sigma_list),
    print "mean:", numpy.mean(sigma_list),
    print "sigma:", numpy.std(sigma_list)

def graph_delta_times(error_list, colour, canvas = None, title=""):
    peak_time_list = [x["x"]*DT for x in error_list] # 0.2 ns per count
    peak_delta_list = [t1-peak_time_list[i] for i, t1 in enumerate(peak_time_list[1:])]
    hist, graph = common.make_root_graph("Peaks", peak_time_list[0:-1], "Time of peak [ns]", peak_delta_list, "Time between peaks [ns]")
    mean_sigma_cut(peak_delta_list)
    if canvas == None:
        canvas = common.make_root_canvas(title+"bpm signal")
        canvas.Draw()
        hist.Draw()
    canvas.cd()
    graph.SetMarkerStyle(6)
    graph.SetMarkerColor(colour)
    graph.Draw("p")
    canvas.Update()    
    return canvas, hist, graph

def graph_peak_errors(error_list, colour, canvas = None, title=""):
    peak_time_list = [x["x"]*DT for x in error_list] # 0.2 ns per count
    peak_error_list = [(x["cov(x,y)"][0][0]**0.5)*DT for x in error_list]
    hist, graph = common.make_root_graph("Peaks", peak_time_list, "Time of peak [ns]", peak_error_list, "Estimated error on peak [ns]")
    if canvas == None:
        canvas = common.make_root_canvas(title+" errors")
        canvas.Draw()
        hist.Draw()
    canvas.cd()
    graph.SetMarkerStyle(6)
    graph.SetMarkerColor(colour)
    graph.Draw("p")
    canvas.Update()    
    return canvas, hist, graph

def graph_peak_times(a_data, errors, colour, y_range, canvas = None, title=""):
    hist, graph = common.make_root_graph("Volts", a_data["time_list"], "Time [ns]", a_data["voltage_list"], "Signal Voltage [V]", ymin=y_range[0], ymax=y_range[1])
    peak_time_list = [x["x"] for x in errors]
    peak_voltage_list = [x["y"] for x in errors]
    peak_hist, peak_graph = common.make_root_graph("Peaks", peak_time_list, "Time of peak [ns]", peak_voltage_list, "Peak voltage")
    if canvas == None:
        canvas = common.make_root_canvas(title+"bpm peaks")
        canvas.Draw()
        hist.Draw()
    canvas.cd()
    peak_graph.SetMarkerColor(colour)
    peak_graph.SetMarkerStyle(4)
    peak_graph.Draw("p")
    graph.SetMarkerColor(colour)
    graph.Draw("p")
    canvas.Update()
    return canvas, hist, graph, peak_hist, peak_graph

def graph_peak_magnitude(pos_peak_with_errors, neg_peak_with_errors, title):
    peak_time_list, peak_voltage_list, peak_err_list = [], [], []
    for i, v_pos in enumerate(pos_peak_with_errors):
        di = 0
        if i >= len(neg_peak_with_errors):
            break
        while neg_peak_with_errors[i+di]["x"] < pos_peak_with_errors[i]["x"]:
            if i+di+1 >= len(neg_peak_with_errors):
                break
            di += 1
        v_neg = neg_peak_with_errors[i+di]
        v_0 = (v_pos["y"]+v_neg["y"])/2.
        v_err = abs(v_0)*(v_pos["cov(x,y)"][1][1]/v_pos["y"]**2+v_neg["cov(x,y)"][1][1]/v_neg["y"]**2)**0.5
        peak_voltage_list.append(v_0)
        peak_err_list.append(v_err*1e3)
        peak_time_list.append(v_pos["x"]*DT)
    print "    Peak voltage mean:", numpy.mean(peak_voltage_list[-20:-10]),
    print "sigma:", numpy.std(peak_voltage_list[-20:-10])
    canvas = common.make_root_canvas(title+" rf voltage")
    canvas.Draw()
    canvas.cd()
    y_range = common.min_max(peak_voltage_list+peak_err_list)
    hist, graph = common.make_root_graph("Peaks", peak_time_list, "Time of peak [ns]", peak_voltage_list, "Peak voltage [kV]", ymin = y_range[0], ymax = y_range[1])
    hist.Draw()
    graph.SetMarkerStyle(6)
    graph.Draw("p")
    hist, graph = common.make_root_graph("Errors", peak_time_list, "Time of peak [ns]", peak_err_list, "Error [V]", ymin = y_range[0], ymax = y_range[1])
    graph.SetMarkerColor(2)
    graph.SetMarkerStyle(6)
    graph.Draw("p")
    canvas.Update()    
    return canvas, hist, graph

def find_peaks(data, window_size):
    print "  Peaks...",
    sys.stdout.flush()
    print 'smoothing...',
    sys.stdout.flush()
    smoothed_data = GaussianSmoothing(window_size/2., window_size, True).smooth(data)
    print 'finding peaks...',
    sys.stdout.flush()
    peak_finder = WindowPeakFinder(window_size, 0., window_size/2)
    peak_list = peak_finder.find_peaks(smoothed_data)
    print 'getting errors for', len(peak_list), 'peaks...',
    sys.stdout.flush()
    error_summary = []
    peak_refiner = RefinePeakFinder(peak_list, 50, 3000, False)
    peak_error_list = peak_refiner.find_peak_errors(data)
    print "found", len(peak_error_list), "errors... Done"
    return peak_error_list

def relative(a_time, time_of_peaks, phase):
    a_time = a_time - phase #time relative to some phase offset
    peak_index = bisect.bisect_left(time_of_peaks, a_time)
    if peak_index > 0:
        return a_time - time_of_peaks[peak_index-1]
    return a_time - time_of_peaks[peak_index]

def find_deltas(data):
    print "Processing Nominal V =", data["RF"]["rf_voltage"]
    data_pair = [data["RF"], data["signal"]]
    delta_t = data_pair[0]["time_list"][100]-data_pair[0]["time_list"][99]
    window_size = int(50/DT)
    times = data_pair[0]["time_list"]
    rf_voltage_list = data_pair[0]["voltage_list"]
    bpm_voltage_list = [-voltage for voltage in data_pair[1]["voltage_list"]]
    rf_peak_to_peak_voltage = max(rf_voltage_list) - min(rf_voltage_list)

    bpm_peak_error_list = find_peaks(bpm_voltage_list, window_size)
    rf_peak_error_list = find_peaks(rf_voltage_list, window_size)
    rf_min_error_list = find_peaks([-x for x in rf_voltage_list], window_size)
    for item in rf_min_error_list:
        item["x"] *= -1.
    data["signal"]["peak_with_errors"] = bpm_peak_error_list
    data["RF"]["pos_peak_with_errors"] = rf_peak_error_list
    data["RF"]["neg_peak_with_errors"] = rf_min_error_list

    bpm_peak_list = [int(x["x"]) for x in bpm_peak_error_list]
    rf_peak_list = [int(x["x"]) for x in rf_peak_error_list]
        
    data["signal"]["peak_indices"] = bpm_peak_list
    data["RF"]["peak_indices"] = rf_peak_list

    print max(bpm_peak_list), len(bpm_voltage_list)
    print max(rf_peak_list), len(rf_voltage_list)
    print "  Calculating deltas"
    bpm_i = 0
    peak_delta_list = []
    peak_time_list = []
    for i, rf_peak in enumerate(rf_peak_error_list[0:-1]):
        rf_time = rf_peak["x"]
        print "RF time", rf_time
        rf_next_time = rf_peak_error_list[i+1]["x"]
        while bpm_i < len(bpm_peak_error_list) and \
              bpm_peak_error_list[bpm_i]["x"] < rf_time:
            bpm_i += 1
        while bpm_i < len(bpm_peak_error_list) and \
              bpm_peak_error_list[bpm_i]["x"] < rf_next_time:
            peak_delta_list.append(bpm_peak_error_list[bpm_i]["x"]-rf_time)
            peak_time_list.append(rf_time)
            bpm_i += 1
    
    data["peak_delta_list"] = [x*DT for x in peak_delta_list]
    data["peak_time_list"] = [x*DT for x in peak_time_list]

    peak_v = numpy.mean([peak["y"] for peak in data["RF"]["pos_peak_with_errors"]])-\
             numpy.mean([peak["y"] for peak in data["RF"]["neg_peak_with_errors"]])
    data["rf_peak_to_peak_voltage"] = peak_v

def plot(plot_dir, data_list_pairs):
    tmp_data = []
    for i, data in enumerate(data_list_pairs):
        pos_peak = numpy.mean([x["y"] for x in data["RF"]["pos_peak_with_errors"]])
        neg_peak = numpy.mean([-x["y"] for x in data["RF"]["neg_peak_with_errors"]])
        voltage = (pos_peak-neg_peak)
        data["rf_peak_to_peak_voltage"] = voltage
        rf_list = ["0.2", "0.32", "0.46", "0.62", "0.81", "3.01", "7.18"]
        if True or str(round(voltage, 2)) in rf_list:
            print "Using", voltage
            tmp_data.append(data)
        else:
            print "Not using", voltage
    data_list_pairs = tmp_data
    peak_delta_canvas = common.make_root_canvas("delta rf vs bpm")
    peak_delta_canvas.Draw()
    peak_delta_hist_canvas = common.make_root_canvas("delta rf vs bpm - hist")
    peak_delta_hist_canvas.Draw()
    graph_list = []
    hist_list = []
    dt_mean_list = []
    dt_err_list = []
    measured_voltage_list = []
    max_y = -1e9
    min_y = 1e9
    delta_hist_min = 1e9
    delta_hist_max = -1e9
    for i, data in enumerate(data_list_pairs):
        """
        voltage_str = "V="+str(round(data["rf_peak_to_peak_voltage"], 2))
        print "Plotting measured "+voltage_str
        formats = ["png", "root"]
        times = data["RF"]["time_list"]

        rf_voltage_list = data["RF"]["voltage_list"]
        rf_peak_error_list = data["RF"]["pos_peak_with_errors"]
        rf_peak_list = [x["x"] for x in rf_peak_error_list]
        rf_period = sum([rf_peak_list[j+1] - rf_peak_list[j] for j, index in enumerate(rf_peak_list[:-1])])  
        rf_period = DT*rf_period/float(len(rf_peak_list))
        rf_frequency = 1./rf_period
        print "  RF Period", rf_period, "Frequency", rf_frequency

        print "  Plotting RF magnitude and error"
        canvas_rf_v, hist_rf_v, graph_rf_v = graph_peak_magnitude(data["RF"]["pos_peak_with_errors"], data["RF"]["neg_peak_with_errors"], title=voltage_str)
        for format in formats:
            canvas_rf_v.Print(plot_dir+voltage_str+"_rf_voltage."+format)

        print "  Plotting peak finding errors"
        canvas_err, hist_err, graph_err = graph_peak_errors(data["signal"]["peak_with_errors"], 4, title=voltage_str)
        graph_peak_errors(data["RF"]["pos_peak_with_errors"], 2, canvas_err)
        for format in formats:
            canvas_err.Print(plot_dir+voltage_str+"_peak_errors."+format)

        print "  Plotting peak to peak distance"
        print "    signal",
        canvas_pp, hist_pp, graph_pp = graph_delta_times(data["signal"]["peak_with_errors"], 4, title=voltage_str)
        print "    RF",
        graph_delta_times(data["RF"]["pos_peak_with_errors"], 2, canvas_pp)
        for format in formats:
            canvas_pp.Print(plot_dir+voltage_str+"_peak_deltas."+format)

        print "  Plotting data"
        y_range = common.min_max(data["RF"]["voltage_list"]+data["signal"]["voltage_list"]+data["dc_signal"]["voltage_list"])
        canvas, hist, graph, peak_hist, peak_graph = graph_peak_times(data["RF"], data["RF"]["pos_peak_with_errors"], 2, y_range, title=voltage_str)
        bpm_peak_indices = data["signal"]["peak_indices"]
        bpm_voltage_list = data["signal"]["voltage_list"]
        bpm_peak_list = [index for index in bpm_peak_indices]
        graph_peak_times(data["signal"], data["signal"]["peak_with_errors"], 4, y_range, canvas)
        graph_peak_times(data["dc_signal"], [{"x":0, "y":0}], 6, y_range, canvas)
        canvas.Update()
        for format in formats:
            canvas.Print(plot_dir+voltage_str+"_signals."+format)
        """

        print "  Plotting distance between bpm and rf peaks"
        peak_delta_canvas.cd()
        peak_time_list = data["peak_time_list"]
        peak_delta_list = data["peak_delta_list"]

        hist, graph = common.make_root_graph(voltage_str, peak_time_list, "time [ns]", peak_delta_list, "dt [ns]")
        min_y = min([hist.GetYaxis().GetXmin(), min_y])
        max_y = max([hist.GetYaxis().GetXmax(), max_y])
        if len(data_list_pairs) > 1:
            hist.Draw()
            color_fraction = float(i)/float(len(data_list_pairs)-1)
        else:
            color_fraction = 1
        color = ROOT.TColor.GetColor(color_fraction, 1.-color_fraction, 0)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.SetMarkerStyle(6)
        graph_list.append(graph)

        peak_delta_hist_canvas.cd()
        filtered_list = []
        for j, peak in enumerate(peak_delta_list):
            if peak_time_list[j] > 30.e3 and \
               peak_time_list[j] < 40.e3:
                filtered_list.append(peak)
        if len(filtered_list) > 0:
            dt_mean_list.append(max(filtered_list))
            #dt_err_list.append(numpy.std(filtered_list)/len(filtered_list)**0.5)
            measured_voltage_list.append(data["rf_peak_to_peak_voltage"])
        hist = common.make_root_histogram(voltage_str, filtered_list, "dt [ns]", 4)
        delta_hist_min = min([hist.GetXaxis().GetXmin(), delta_hist_min])
        delta_hist_max = max([hist.GetXaxis().GetXmax(), delta_hist_max])
        hist.SetLineColor(color)
        hist_list.append(hist)

    peak_delta_canvas.cd()
    hist, graph = common.make_root_graph(voltage_str, peak_time_list, "t [ns]", peak_delta_list, "dt [ns]", ymin=0., ymax=700.)
    hist.Draw() # redraw hist to get axes right
    for graph in graph_list:
        graph.Draw("p")
    legend = common.make_root_legend(peak_delta_canvas, graph_list)
    legend.Draw()
    peak_delta_canvas.Update() 
    for format in formats:
        peak_delta_canvas.Print(plot_dir+"bpm_to_rf_deltas."+format)

    peak_delta_hist_canvas.cd()
    print delta_hist_min, delta_hist_max
    hist = common.make_root_histogram(voltage_str, [-1], "dt [ns]", 1000, [-1], "Counts", 1000, xmin=delta_hist_min, xmax=delta_hist_max, ymin=1e-2, ymax=20.)
    hist.Draw() # redraw hist to get axes right
    for a_hist in hist_list:
        a_hist.Draw("same")
    legend = common.make_root_legend(peak_delta_hist_canvas, hist_list)
    legend.Draw()
    peak_delta_hist_canvas.Update()
    for format in formats:
        peak_delta_hist_canvas.Print(plot_dir+"bpm_to_rf_deltas_hist."+format)

    if len(dt_mean_list) < 2:
        return
    phase_mean_list = [2.*math.pi*dt/rf_period for dt in dt_mean_list[1:]]
    measured_voltage_list = measured_voltage_list[1:]
    for popper in []:
        phase_mean_list.pop(popper)
        phase_err_list.pop(popper)
        measured_voltage_list.pop(popper)
    synchronous_phase_canvas = common.make_root_canvas("Synchronous phase vs voltage")
    hist, graph = common.make_root_graph("Synchronous phase vs voltage", phase_mean_list, "#phi [rad]", measured_voltage_list, "V [au]")
    graph = ROOT.TGraph(len(phase_mean_list))
    graph.SetMarkerStyle(4)
    for i in range(len(phase_mean_list)):
        graph.SetPoint(i, phase_mean_list[i], measured_voltage_list[i])

    print phase_mean_list 
    print measured_voltage_list
    hist.Draw()
    graph.Draw('p')
    print "Doing unweighted fit"
    fit_unweighted = ROOT.TF1("fit", "[0]/sin([1]+x)")
    fit_unweighted.SetParameter(0, 0.1)
    fit_unweighted.SetParameter(1, 0.)
    fit_unweighted.SetLineColor(4)
    graph.Fit(fit_unweighted, "EX0")
    fit_unweighted.Draw("SAME")
    synchronous_phase_canvas.Update()    
    for format in formats:
        synchronous_phase_canvas.Print(plot_dir+"synchronous_phase_vs_voltage_unweighted."+format)

def main(args):
    if '--plot' not in args and '--process' not in args:
        raise RuntimeError("Should specify --process and/or --plot on command line")

    if '--process' in args:
        try:
            print "Using data from", os.environ["run_2014_06"]
        except KeyError:
            print "\nMissing environment variable"
            raise
        file_index_string = open('file_index.json').read()
        file_index = json.loads(file_index_string)
        filename = "data_output"
        try:
            index_number = int(args[args.index('--process')+1])
        except IndexError:
            for index in range(len(file_index)):
                job = ['bsub', '-q', 'scarf-ibis', '-W', '02:00', 'python', 'phi_s.py', '--process', str(index)]
                print ' '.join(job)
                proc = subprocess.Popen(job)
            return
        fout = open(filename+"_"+str(index_number)+'.json', 'w')
        data = load_data_from_listing(file_index[index_number])
        find_deltas(data)
        print >> fout, json.dumps(data)

    if '--plot' in args:
        filename = "data_output_all.json"
        print "Reading from", filename
        fin = open(filename)
        data_list_pairs = []
        lines = [line for i, line in enumerate(fin.readlines())]
        for line in lines:
            try:
                data_list_pairs.append(json.loads(line))
            except ValueError:
                sys.excepthook(*sys.exc_info())
        plot("plots/", data_list_pairs)
        raw_input()


if __name__ == "__main__":
    main(sys.argv)

