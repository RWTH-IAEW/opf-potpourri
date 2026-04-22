"""AC-OPF analysis for representative simbench load cases.

Workflow
--------
This script demonstrates two successive AC-OPF steps on the simbench LV
rural1 network:

**Step 1 – reactive power optimisation (single snapshot)**

1. Load the network and select a specific time-step profile
   (``select_profile_idx = 1190``).
2. Derive reactive power limits for all PV generators from a 0.95
   power-factor envelope.
3. Fix active power dispatch and run a reference pandapower power flow.
4. Constrain the external grid active power to the pandapower result (so the
   OPF only dispatches reactive power from the generators).
5. Solve an AC-OPF with ``add_reactive_power_flow_objective`` – minimises the
   sum of squared reactive power injections from static generators.
6. Print the generator set-points before and after the optimisation.

**Step 2 – voltage-deviation minimisation per load case**

For each standard simbench load case (``"lW"`` – low wind, ``"hL"`` – high
load):

1. Scale loads and renewable generation according to the case factors.
2. Set the slack bus voltage to the case-specific value.
3. Solve an AC-OPF with ``add_voltage_deviation_objective`` – minimises the
   sum of squared per-unit voltage deviations from 1.0.
4. Print the ext-grid active / reactive power dispatch and the active power on
   line 0.

Requires IPOPT (installed via conda-forge or environment.yaml).
"""

import copy
import warnings

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # Network setup                                                        #
    # ------------------------------------------------------------------ #
    net = sb.get_simbench_net("1-LV-rural1--0-sw")

    # Select a single time-step snapshot from the simbench time-series.
    profiles = sb.get_absolute_values(
        net, profiles_instead_of_study_cases=True
    )
    select_profile_idx = 1190

    net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[select_profile_idx]
    net.load["p_mw"] = profiles[("load", "p_mw")].iloc[select_profile_idx]
    net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[select_profile_idx]

    # Reactive power envelope for PV: Q-limit derived from 0.95 power factor.
    power_factor = 0.95
    q_max = np.sqrt(
        (net.sgen["p_mw"] / power_factor) ** 2 - net.sgen["p_mw"] ** 2
    )
    net.sgen["max_q_mvar"] = q_max
    net.sgen["min_q_mvar"] = -q_max

    # Fix active power at current dispatch; mark generators as controllable
    # so the OPF can optimise their reactive power.
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = net.sgen["p_mw"]
    net.sgen["controllable"] = True

    # Reference power flow to obtain current ext-grid operating point.
    pp.runpp(net, voltage_depend_loads=False)

    # ------------------------------------------------------------------ #
    # Step 1: reactive-power minimisation (single snapshot)               #
    # ------------------------------------------------------------------ #
    net_case_opf = copy.deepcopy(net)

    # Pin the external grid to the pandapower operating point; allow only a
    # small reactive-power band so the OPF must draw Q from the generators.
    net_case_opf.ext_grid["max_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
    net_case_opf.ext_grid["min_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
    net_case_opf.ext_grid["max_q_mvar"] = 0.05
    net_case_opf.ext_grid["min_q_mvar"] = -0.05

    hc = ACOPF(net_case_opf)
    hc.add_OPF()
    hc.add_reactive_power_flow_objective()

    print("Grid-state before OPF (from net.sgen):")
    print(net.sgen[["p_mw", "q_mvar"]])

    print("Pyomo initial set-points (PsG param, QsG param):")
    print("  SGEN P [pu]:")
    for g in hc.model.sG:
        print(f"    sgen {g}: {pyo.value(hc.model.PsG[g]):.4f}")
    print("  SGEN Q [pu]:")
    for g in hc.model.sG:
        print(f"    sgen {g}: {pyo.value(hc.model.QsG[g]):.4f}")

    # hc.solve(solver="neos", print_solver_output=True)
    hc.solve(solver="ipopt", print_solver_output=True)

    print("Pyomo optimised set-points (psG var, qsG var):")
    print("  SGEN P [pu]:")
    for g in hc.model.sG:
        print(f"    sgen {g}: {pyo.value(hc.model.psG[g]):.4f}")
    print("  SGEN Q [pu]:")
    for g in hc.model.sG:
        print(f"    sgen {g}: {pyo.value(hc.model.qsG[g]):.4f}")

    # ------------------------------------------------------------------ #
    # Step 2: voltage-deviation minimisation per load case                #
    # ------------------------------------------------------------------ #
    case_keys = ["lW", "hL"]

    hcs = []  # solved ACOPF objects, one per load case
    obj = []  # objective values per load case

    for case in case_keys:
        net_case = copy.deepcopy(net)

        factors = net_case.loadcases.loc[case]

        # Scale loads and renewables according to the load-case factors.
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

        # No add_OPF() here – the intent is to find the reactive power dispatch
        # that minimises voltage deviation over the free AC power flow, without
        # imposing generator capacity or line loading limits.
        hc = ACOPF(net_case)
        hc.add_voltage_deviation_objective()
        # hc.solve(solver="neos")
        hc.solve(solver="ipopt")

        hcs.append(copy.deepcopy(hc))
        obj.append(pyo.value(hc.model.obj_v_deviation))

        print(f"\nLoad case '{case}' – ext-grid dispatch [pu]:")
        for g in hc.model.eG:
            print(
                f"  ext_grid {g}  P = {pyo.value(hc.model.pG[g]):.4f}"
                f"  Q = {pyo.value(hc.model.qG[g]):.4f}"
            )

        line_flow = pyo.value(hc.model.pLfrom[0])
        print(f"  Line 0 active power (from-end) [pu]: {line_flow:.4f}")
