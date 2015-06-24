import bisect
import math
import json
import sys

import numpy
import ROOT

import xboa.common as common
import xboa.common.root_wrapper
import xboa.algorithms.peak_finder
import xboa.algorithms.smoothing

def periodic(a_time, period, phase):
    a_time = a_time - phase #time relative to first RF peak
    base_time = math.floor(a_time/period)*period
    a_time = a_time - base_time
    return a_time

def relative(a_time, time_of_peaks, phase):
    a_time = a_time - phase #time relative to some phase offset
    peak_index = bisect.bisect_left(time_of_peaks, a_time)
    if peak_index > 0:
        return a_time - time_of_peaks[peak_index-1]
    return a_time - time_of_peaks[peak_index]

def get_period(data):
    rf_peak_indices = data["RF"]["peak_indices"]
    rf_peak_list = [data["RF"]["time_list"][i]/1e-9 for i in rf_peak_indices]
    rf_peak_deltas_list = [None]*(len(rf_peak_list)-1)

    for i, rf_time in enumerate(rf_peak_list[0:-1]):
        rf_next_time = rf_peak_list[i+1]
        rf_peak_deltas_list[i] = rf_next_time-rf_time
    canvas = common.make_root_canvas("rf_period")
    canvas.Draw()
    hist = common.make_root_histogram("times", rf_peak_deltas_list, "RF period [ns]", 100)
    hist.Draw()
    hist.SetStats(True)
    canvas.Update()
    canvas.Draw()

    return numpy.mean(rf_peak_deltas_list), rf_peak_list[0], canvas

def do_fit(graph, xmin, xmax, ymin, ymax):
    graph_fitted = ROOT.TGraph()
    fit_parameters = {"fit_parameters":[]}
    x_list, y_list = [], []
    for i in range(graph.GetN()):
        x = ROOT.Double()
        y = ROOT.Double()
        graph.GetPoint(i, x, y)
        if x > xmin and x < xmax and y > ymin and y < ymax:
            x_list.append(x)
            y_list.append(y)
    if len(x_list) > 5:
        fit_parameters["mean_dt"] = numpy.mean(y_list)
        fit_parameters["std_dt"] = numpy.std(y_list)
        fit_parameters["n_dt"] = len(y_list)
        fit_parameters["delta_dt"] = (max(y_list)-min(y_list))/2.
        fit_parameters["peak_dt"], peaks_canvas = peak_hist(y_list)
    else:
        fit_parameters["mean_dt"] = None
        fit_parameters["std_dt"] = None
        fit_parameters["peak_dt"] = None
    fit_parameters["peak_dt_min"] = min([a_peak["lower_bound"] for a_peak in fit_parameters["peak_dt"]])
    fit_parameters["peak_dt_max"] = max([a_peak["upper_bound"] for a_peak in fit_parameters["peak_dt"]])
    n_points = 0
    n_points_fitted = 0
    for i in range(graph.GetN()):
        x = ROOT.Double()
        y = ROOT.Double()
        graph.GetPoint(i, x, y)
        if x > xmin and x < xmax:
            n_points += 1
            if fit_parameters["peak_dt_min"] < fit_parameters["peak_dt_max"]:
                if y > fit_parameters["peak_dt_min"] and y < fit_parameters["peak_dt_max"]:
                    n_points_fitted += 1
                    graph_fitted.SetPoint(i, x, y)
            else: # wrapped around...
                if y > fit_parameters["peak_dt_min"] or y < fit_parameters["peak_dt_max"]:
                    n_points_fitted += 1
                    graph_fitted.SetPoint(i, x, y)
    formula = "[0]*sin(6.2831853*x/[1]+[2])+[3]"
    fit = ROOT.TF1("fit", formula)
    fit.SetParameter(0, 50)
    fit.SetParameter(1, -200.e3)
    fit.SetParameter(2, -150.e3)
    fit.SetParameter(3, 350.)
    fit.SetLineColor(4)  
    #graph_fitted.Fit(fit, "", "", xmin, xmax)
    graph_fitted.SetMarkerStyle(6)
    graph_fitted.SetMarkerColor(7)
    xboa.common.root_wrapper.keep_root_object(graph_fitted)
    xboa.common.root_wrapper.keep_root_object(fit)
    for i in range(4):
        fit_parameters["fit_parameters"].append([fit.GetParameter(i), fit.GetParError(i)])
    fit_parameters["chi2"] = fit.GetChisquare()
    fit_parameters["n_dof"] = fit.GetNDF()
    fit_parameters["fit_function"] = formula
    fit_parameters["cut"] = {"xmin":xmin, "xmax":xmax, "ymin":ymin, "ymax":ymax}
    fit_parameters["fitted_fraction"] = float(n_points_fitted)/float(n_points)

    return graph_fitted, fit, peaks_canvas, fit_parameters

def peak_hist(peak_list):
    canvas = common.make_root_canvas("projected peaks")
    hist = common.make_root_histogram("projected peaks", peak_list, "#deltat [ns]", 200, xmin=0, xmax=630.)
    bins = []
    for bin_index in range(1, hist.GetNbinsX()+1):
        bins.append(hist.GetBinContent(bin_index))
    smoothing = xboa.algorithms.smoothing.GaussianSmoothing(2., 3, True)
    smoothed = smoothing.smooth(bins)
    finder = xboa.algorithms.peak_finder.WindowPeakFinder(10, -1, 1)
    peak_graph_indices, peak_x, peak_y = [], [], []
    peaks = []
    for peak in finder.find_peaks(smoothed):
        a_peak = {"peak":hist.GetBinCenter(peak+1), "magnitude":bins[peak], "err":0.}
        if a_peak["magnitude"] < 0.5*max(smoothed):
            continue
        peak_graph_indices.append(peak)
        for bin in range(peak, 0, -1)+range(len(bins)-1, peak, -1):
            if smoothed[bin] < smoothed[peak]/2.:
                peak_graph_indices.append(bin)
                a_peak["lower_bound"] = hist.GetBinCenter(bin+1)
                break
        for bin in range(peak, len(bins), +1)+range(0, peak, 1):
            if smoothed[bin] < smoothed[peak]/2.:
                peak_graph_indices.append(bin)
                a_peak["upper_bound"] = hist.GetBinCenter(bin+1)
                break
        peaks.append(a_peak)
    for peak_graph_index in peak_graph_indices:
        peak_x.append(hist.GetBinCenter(peak_graph_index+1))
        peak_y.append(smoothed[peak_graph_index])
    dummy, graph = xboa.common.make_root_graph("peaks", peak_x, "", peak_y, "")
    graph.SetMarkerStyle(24)
    graph.SetMarkerColor(2)
    canvas.Draw()
    hist.Draw()
    graph.Draw('p')
    canvas.Update()
    return peaks, canvas


def normalise_volts(volts, periodic_times):
    test_time = -1e9
    norm_volts = []
    test_volts = []
    for i, a_volt in enumerate(volts):
        if periodic_times[i] < test_time and len(test_volts) > 1:
            min_test_volt = min(test_volts)
            max_test_volt = max(test_volts)
            delta_volts = (max_test_volt-min_test_volt)
            test_volts = [(a_volt_2-min_test_volt)/delta_volts for a_volt_2 in test_volts]
            norm_volts += test_volts
            test_volts = []
        test_time = periodic_times[i]
        test_volts.append(a_volt)
    norm_volts += test_volts
            
    #print "PERIOD", periodic_times
    #print "VOLTS", test_volts
    #print norm_volts[-3000:]
    print "."
    return norm_volts

def plot_one(plot_dir, data, fout): 
    print "    Parsing data"
    time_of_peaks = sorted([x["x"]*0.2 for x in data["RF"]["pos_peak_with_errors"]])

    print "  Plotting distance between bpm and rf peaks"
    peak_time_list = data["peak_time_list"]
    peak_delta_list = data["peak_delta_list"]

    dummy, graph = common.make_root_graph("peaks", peak_time_list, "time [ns]", peak_delta_list, "dt [ns]")
    graph.SetMarkerStyle(7)

    graph_fitted, fit, peaks_canvas, fit_parameters = do_fit(graph, 150.e3, 400.e3, 0., 700.)

    pos_peak = numpy.mean([x["y"] for x in data["RF"]["pos_peak_with_errors"]])
    neg_peak = numpy.mean([-x["y"] for x in data["RF"]["neg_peak_with_errors"]])
    peak_to_peak_voltage = (pos_peak-neg_peak)

    volts = [-a_volt for a_volt in data["signal"]["voltage_list"][::]]
    times = [a_time/1e-9 for a_time in data["signal"]["time_list"]]
    time_zero = times[0]
    times = [a_time - time_zero for a_time in times]

    print "    Getting period"
    period, phase, frequency_canvas = get_period(data)

    print "    Getting periodic times"
    periodic_times = [relative(a_time, time_of_peaks, times[0]) for a_time in times]
    print "    Voltage:", peak_to_peak_voltage, "T:", period, "N_events:", len(periodic_times)
    norm_volts = normalise_volts(volts, periodic_times)
    canvas = common.make_root_canvas("times")
    min_times, max_times = min(time_of_peaks), max(time_of_peaks)
    n_bins = len(time_of_peaks)-1
    print "    NBins", int(n_bins), "time min/max", min_times, max_times
    hist = ROOT.TH2D("", ";Time in cycle [ns*1e3];#deltat [ns]", int(n_bins), min_times, max_times, int(period), 0., int(period))
    for i, a_time in enumerate(times):
        bin = hist.FindFixBin(a_time, periodic_times[i])
        if abs(hist.GetBinContent(bin)) < 1e-9:
            hist.SetBinContent(bin, norm_volts[i])
    hist.SetStats(False)
    hist.SetTitle("Voltage: "+str(round(peak_to_peak_voltage/2., 2)))
    xboa.common.root_wrapper.keep_root_object(hist)
    canvas.Draw()
    #canvas.SetLogz()
    hist.Draw("COLZ")
    graph.Draw('p')
    graph_fitted.Draw('p')

    canvas.Update()


    name = plot_dir+"/V="+str(round(peak_to_peak_voltage, 2))
    name = name.replace(".", "_")
    for format in ["png", "root"]:
        canvas.Print(name+"_fitted_bpm_to_rf_deltas."+format)
        frequency_canvas.Print(name+"_frequency_distribution."+format)
        peaks_canvas.Print(name+"_peaks."+format)
    fit_parameters["peak_to_peak_voltage"] = peak_to_peak_voltage
    fit_parameters["rf_period"] = period
    print "    Fit parameters", json.dumps(fit_parameters, indent=2)

    print >> fout, json.dumps(fit_parameters)

def main():
    fin = open("data_output_ref.json")
    fout = open("data_summary.json", "w")
    for i, line in enumerate(fin.readlines()):
        print "Loading next line...", i
        ROOT.gROOT.SetBatch(i > 2)
        try:
            data = json.loads(line)
            plot_one("plots/", data, fout)
        except ValueError:
            sys.excepthook(*sys.exc_info())

if __name__ == "__main__":
    main()
    print "Finished - press <CR>"
    raw_input()

