import numpy as np
import simbench as sb
import pyomo.environ as pyo
import pandapower as pp

from src.potpourri.models.ACOPF_base import ACOPF

if __name__ == "__main__":
    # net = sb.get_simbench_net("1-HV-mixed--0-no_sw")
    net = sb.get_simbench_net("1-LV-rural1--0-sw")

    case_keys = ["lW", "hL"]

    hcs = []
    obj = []

    # Load net data
    profiles = sb.get_absolute_values(
        net, profiles_instead_of_study_cases=True
    )

    select_profile_idx = 1190

    net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[select_profile_idx]
    net.load["p_mw"] = profiles[("load", "p_mw")].iloc[select_profile_idx]
    net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[select_profile_idx]

    power_factor = 0.95
    q_max = np.sqrt(
        (net.sgen["p_mw"] / power_factor) ** 2 - net.sgen["p_mw"] ** 2
    )
    net.sgen["max_q_mvar"] = q_max
    net.sgen["min_q_mvar"] = -q_max

    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = net.sgen["p_mw"]
    net.sgen["controllable"] = True
    pp.runpp(net)

    net_case_opf = net.deepcopy()

    net_case_opf.ext_grid["max_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
    net_case_opf.ext_grid["min_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
    net_case_opf.ext_grid["max_q_mvar"] = 0.05
    net_case_opf.ext_grid["min_q_mvar"] = -0.05

    hc = ACOPF(net_case_opf)
    hc.add_OPF()
    hc.add_reactive_power_flow_objective()
    print("Grid-state before OPF:")
    print(net.sgen[["p_mw", "q_mvar"]])
    print("Grid-state before OPF:")
    print("SGEN P:")
    for g in hc.model.sG:
        print(pyo.value(hc.model.PsG[g]))
    print("SGEN Q:")
    for g in hc.model.sG:
        print(pyo.value(hc.model.QsG[g]))

    # hc.solve(solver='neos', print_solver_output=True)
    hc.solve(solver="ipopt", print_solver_output=False)

    # Print summary of changes
    print("Grid-state after OPF:")
    print(net.sgen[["p_mw", "q_mvar"]])
    print("Grid-state after OPF:")
    print("SGEN P:")
    for g in hc.model.sG:
        print(pyo.value(hc.model.psG[g]))
    print("SGEN Q:")
    for g in hc.model.sG:
        print(pyo.value(hc.model.qsG[g]))

    for case in case_keys:
        net_case = net.deepcopy()

        factors = net_case.loadcases.loc[case]

        net_case.load.p_mw *= factors["pload"]
        net_case.load.q_mvar *= factors["qload"]

        net_case.sgen.loc[net_case.sgen.type == "Wind", "scaling"] = factors[
            "Wind_p"
        ]
        net_case.sgen.loc[net_case.sgen.type == "PV", "scaling"] = factors[
            "PV_p"
        ]
        net_case.sgen.loc[
            (net_case.sgen.type != "Wind") & (net_case.sgen.type != "PV"),
            "scaling",
        ] = factors["RES_p"]

        net_case.ext_grid.vm_pu = factors["Slack_vm"]

        hc = ACOPF(net_case)
        hc.add_voltage_deviation_objective()
        hc.solve(solver="neos")

        # hcs.append(copy.deepcopy(hc))
        # obj.append(pyo.value(hc.model.obj))

        # for g in hc.model.sG:
        #    print(f"Gen {g}: {pyo.value(hc.model.pG[g])}")
        #    print(f"Gen {g}: {pyo.value(hc.model.qG[g])}")

        line_flow = pyo.value(hc.model.pLfrom[0])
        print(f"Line flow: {line_flow}")
