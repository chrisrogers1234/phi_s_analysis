import bisect
import math
import json

import numpy
import ROOT

import xboa.common
import xboa.common.root_wrapper


def periodic(a_time, period, phase):
    a_time = a_time - phase #time relative to first RF peak
    base_time = math.floor(a_time/period)*period
    a_time = a_time - base_time
    return a_time

def get_parameters():
    global MINUIT
    optimiser = MINUIT
    phi_s = ROOT.Double(0.) # numpy.array((1,), dtype=numpy.float64)
    v_0 = ROOT.Double(0.) # numpy.array((1,), dtype=numpy.float64)
    error = ROOT.Double(0.) # numpy.array((1,), dtype=numpy.float64)
    optimiser.GetParameter(0, v_0, error)
    optimiser.GetParameter(1, phi_s, error)
    return v_0, phi_s

def get_chi2(v_0, phi_s):
    global DELTA_PHI_LIST, VOLTS_LIST, DELTA_PHI_ERR_LIST, VOLTS_ERR_LIST
    phi0 = -phi_s+1e-9 #min([dphi - DELTA_PHI_ERR_LIST[i] for i, dphi in enumerate(DELTA_PHI_LIST)])
    phi1 = math.pi/2. - phi_s #max([dphi + DELTA_PHI_ERR_LIST[i] for i, dphi in enumerate(DELTA_PHI_LIST)])
    test_phi_list = [phi_test*(phi1-phi0)/1000.+phi0 for phi_test in range(1000, 0, -1)]
    test_volts_list = [v_0/math.sin(phi_test+phi_s) for phi_test in test_phi_list]
    test_list = sorted(zip(test_volts_list, test_phi_list))
    graph = ROOT.TGraph(1000)
    for i in range(1000):
        graph.SetPoint(i, test_list[i][1], test_list[i][0])
    graph.SetMarkerColor(4)
    graph.SetMarkerStyle(6)
    data = {"runs":[], "score":0., "phi0":phi0, "phi1":phi1, "graph":graph}
    for i, ref_volts in enumerate(VOLTS_LIST):
        ref_phi = DELTA_PHI_LIST[i]
        ref_volts_err = VOLTS_ERR_LIST[i]
        ref_phi_err = DELTA_PHI_ERR_LIST[i]
        chi2_list = []
        for test_volts, test_phi in test_list:
            miss_volts = (test_volts-ref_volts)/ref_volts_err
            miss_phi = (test_phi-ref_phi)/ref_phi_err
            chi2_list.append(miss_volts**2 + miss_phi**2)
        data["score"] += min(chi2_list)
    return data

def minuit_function(n_pars, pars, score, buff, err):
    v_0, phi_s = get_parameters()
    score[0] = get_chi2(v_0, phi_s)["score"]
    #print v_0, phi_s, "**", score[0]

def do_minuit(delta_phi_list, volts_list, delta_phi_err_list, volts_err_list, volts_min):
    global MINUIT, DELTA_PHI_LIST, VOLTS_LIST, DELTA_PHI_ERR_LIST, VOLTS_ERR_LIST
    optimiser = ROOT.TMinuit()
    optimiser.SetPrintLevel(-1)
    optimiser.DefineParameter(0, "V0", min(volts_list), 0.1, 0.01, max(volts_list))
    optimiser.DefineParameter(1, "phi_s", -3.5, 0.1, -5., -1.)
    optimiser.SetFCN(minuit_function)
    MINUIT = optimiser
    DELTA_PHI_LIST = delta_phi_list
    VOLTS_LIST = volts_list
    DELTA_PHI_ERR_LIST = delta_phi_err_list
    VOLTS_ERR_LIST = volts_err_list
    args = numpy.array([500, 1e-6], dtype=numpy.float64)
    ierr = ROOT.Long()
    optimiser.mnexcm("SIMPLEX", args, 2, ierr)
    v_0, phi_s = get_parameters()
    return v_0, phi_s

def do_fit(volts_list, delta_phi_list, delta_phi_err_list, color, volts_min, volts_max):
    graph_fit = ROOT.TGraphErrors()
    graph_plot = ROOT.TGraphErrors()
    volts_min_index = -1
    for i, volts in enumerate(volts_list):
        graph_plot.SetPoint(i, delta_phi_list[i], volts)
        graph_plot.SetPointError(i, delta_phi_err_list[i], 0.005)
        if volts < volts_min or volts > volts_max:
            volts_min_index = i
            continue
        graph_fit.SetPoint(i, delta_phi_list[i], volts)
        graph_fit.SetPointError(i, delta_phi_err_list[i], 0.005)
        #print "    V:", volts, "dphi:", delta_phi_list[i], "err(dphi):", delta_phi_err_list[i]
    graph_plot.SetMarkerColor(color)
    graph_plot.SetMarkerStyle(25)
    graph_plot.Draw('p')
    graph_fit.Draw('p')
    v_0, phi_s = do_minuit(delta_phi_list[volts_min_index:], volts_list[volts_min_index:], delta_phi_err_list[volts_min_index:], [0.005 for i in volts_list[volts_min_index:]], volts_list[volts_min_index-2])
    ndof = len(delta_phi_list[volts_min_index:])
    chi2 = get_chi2(v_0, phi_s)
    fit = get_chi2(v_0, phi_s)["graph"]
    fit.Draw('p')
    graph_plot.Draw('p')
    xboa.common.root_wrapper.keep_root_object(fit)
    xboa.common.root_wrapper.keep_root_object(graph_fit)
    xboa.common.root_wrapper.keep_root_object(graph_plot)
    print "V_min", round(volts_list[volts_min_index], 5), "V_0", round(v_0, 5), "phi_s", round(phi_s, 5), "NDof", ndof,  "Chi2", round(chi2['score'], 5), "Prob", round(ROOT.TMath.Prob(chi2["score"], ndof), 5)
    return (volts_list[volts_min_index], v_0, phi_s, ndof,  chi2['score'], ROOT.TMath.Prob(chi2["score"], ndof))

def plot_phi_s(data):

    print "Fit 1"
    delta_phi, delta_phi_err, volts, fitted_fractions = [], [], [], []
    for i, a_datum in enumerate(data):
        fitted_fractions.append(a_datum["fitted_fraction"])
        volts.append(a_datum["peak_to_peak_voltage"]/2.)
        rf_period = a_datum["rf_period"]
        dt = (a_datum["peak_dt_max"]+a_datum["peak_dt_min"])/2.
        if a_datum["peak_dt_max"] > a_datum["peak_dt_min"]:
            dt_err = (a_datum["peak_dt_max"]-a_datum["peak_dt_min"])/2.
        else:
            dt_err = (a_datum["peak_dt_max"]+(632.-a_datum["peak_dt_min"]))/2.
        delta_phi.append(math.pi*2.*dt/rf_period)
        delta_phi_err.append(math.pi*2.*dt_err/rf_period)
    make_fitted_fractions_graph(volts, fitted_fractions)
    fit_values = []
    canvas = xboa.common.make_root_canvas("phi_s")
    hist, dummy = xboa.common.make_root_graph("", [0.], "Bunch phase [rad]", [0.], "RF Voltage [V]", xmin=0., xmax=math.pi*2., ymin=0., ymax=8.)
    hist.Draw()
    for v_min in volts[2:12]:
        fit_values.append(do_fit(volts, delta_phi, delta_phi_err, 2, v_min-1e-6, 8.0))
    print "Lowest voltage setting included [kV] & Fitted energy loss [keV] & Fitted phase [rad] & Number of degrees of freedom & $\\sum \\chi^2$ & Probability & Fraction of hits in cuts\\\\"
    for fit in fit_values:
        a_fitted_fraction = fitted_fractions[volts.index(fit[0])]
        for value in fit[:-1]:
            print round(value, 5), "&",
        print round(fit[-1], 2), "&", a_fitted_fraction, "\\\\"
    canvas.Update()
    canvas.Print("plots/synchronous_phase_vs_voltage.root")
    canvas.Print("plots/synchronous_phase_vs_voltage.png")

def make_fitted_fractions_graph(volts, fitted_fractions):
    canvas = xboa.common.make_root_canvas("Fitted Fractions", "Fitted Fractions")
    for x, y in zip(volts, fitted_fractions):
        print x, y
    hist, graph = xboa.common.make_root_graph("Fitted Fractions", fitted_fractions, "Fitted fraction", volts, "RF Voltage [kV]")
    canvas.Draw()
    hist.Draw()
    graph.Draw()
    canvas.Update()


def main():
    fin = open("data_summary.ref")
    data = sorted([json.loads(line) for line in fin.readlines()], key = lambda x: x["peak_to_peak_voltage"])
    plot_phi_s(data)
    

if __name__ == "__main__":
    main()
    print "Finished - press <CR> to continue"
    raw_input()

