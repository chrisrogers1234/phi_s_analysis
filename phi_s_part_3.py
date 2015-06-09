import math
import json

import ROOT

import xboa.common
import xboa.common.root_wrapper


def periodic(a_time, period, phase):
    a_time = a_time - phase #time relative to first RF peak
    base_time = math.floor(a_time/period)*period
    a_time = a_time - base_time
    return a_time

def print_voltages(data):
    volts = [a_datum['peak_to_peak_voltage']/2. for a_datum in data]
    for a_volt in sorted(volts):
        print "\item", round(a_volt, 3)

def do_fit(volts_list, delta_phi_list, delta_phi_err_list, color, volts_min, volts_max):
    graph_fit = ROOT.TGraphErrors()
    graph_plot = ROOT.TGraphErrors()
    for i, volts in enumerate(volts_list):
        graph_plot.SetPoint(i, delta_phi_list[i], volts)
        graph_plot.SetPointError(i, delta_phi_err_list[i], 0.005)
        if volts < volts_min or volts > volts_max:
            continue
        graph_fit.SetPoint(i, delta_phi_list[i], volts)
        graph_fit.SetPointError(i, delta_phi_err_list[i], 0.005)
        print "    V:", volts, "fit dphi:", delta_phi_list[i], "err(dphi):", delta_phi_err_list[i]
    graph_plot.SetMarkerColor(color)
    graph_plot.SetMarkerStyle(25)
    graph_plot.Draw('p')
    graph_fit.Draw('p')
    fit = ROOT.TF1("fit", "[0]/sin(x+[1])")
    fit.SetLineColor(color)
    fit.SetParLimits(0, 0.835302652144/2., 1.02319362021/2.)
    fit.SetParameter(1, 3.0)
    graph_fit.Fit(fit, "EX0")
    fit.Draw('SAME')
    xboa.common.root_wrapper.keep_root_object(graph_fit)
    xboa.common.root_wrapper.keep_root_object(graph_plot)
    xboa.common.root_wrapper.keep_root_object(fit)
    

def plot_phi_s(data):
    canvas = xboa.common.make_root_canvas("phi_s")
    hist, dummy = xboa.common.make_root_graph("", [0.], "Bunch phase [rad]", [0.], "RF Voltage [V]", xmin=2., xmax=6.0, ymin=0., ymax=4.)
    hist.Draw()

    print "Fit 1"
    delta_phi, delta_phi_err, volts = [], [], []
    for i, a_datum in enumerate(data):
        volts.append(a_datum["peak_to_peak_voltage"]/2.)
        rf_period = a_datum["rf_period"]
        delta_phi.append(math.pi*2.*a_datum["fit_parameters"][3][0]/rf_period)
        delta_phi_err.append(math.pi*2.*a_datum["fit_parameters"][3][1]/rf_period)
#    do_fit(volts, delta_phi, delta_phi_err, 4, 1.1, 4.0)
      
    print "Fit 2"
    delta_phi, delta_phi_err = [], []
    for i, a_datum in enumerate(data):
        rf_period = a_datum["rf_period"]
        delta_phi.append(math.pi*2.*a_datum["mean_dt"]/rf_period)
        delta_phi_err.append(math.pi*2.*a_datum["std_dt"]/rf_period)
    do_fit(volts, delta_phi, delta_phi_err, 2, 0.4, 4.0)

    delta_phi = []
    for i, a_datum in enumerate(data):
        rf_period = a_datum["rf_period"]
        magnitude = -1.
        #print i, volts[i], ":",
        for dat in a_datum["peak_dt"]:
            #print dat['magnitude'], math.pi*2.*dat['peak']/rf_period, "**",
            if dat['magnitude'] > magnitude:
                magnitude = dat['magnitude']
                max_dphi = math.pi*2.*dat['peak']/rf_period
        #print "Done", magnitude, max_dphi
        delta_phi.append(max_dphi)
    delta_phi[17] = 3.79013655584
#    do_fit(volts, delta_phi, delta_phi_err, 8, 0.5, 2.5)

    canvas.Update()
    canvas.Print("plots/synchronous_phase_vs_voltage.root")
    canvas.Print("plots/synchronous_phase_vs_voltage.png")

def main():
    fin = open("data_summary.json")
    data = sorted([json.loads(line) for line in fin.readlines()], key = lambda x: x["peak_to_peak_voltage"])
    plot_phi_s(data)
    

if __name__ == "__main__":
    main()
    print "Finished - press <CR> to continue"
    raw_input()

