import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np


class Basemodel:
    def __init__(self, net):

        self.net = copy.deepcopy(net)
        self.net.gen.index += len(self.net.sgen.index)

        self.gen_set = self.net.gen.index
        self.sgen_set = self.net.sgen.index
        self.generators = pd.concat([self.net.sgen, self.net.gen])
        self.gen_all_set = self.generators.index

        # --- Set indices ---
        self.bus_set = net.bus.index
        self.demand_set = net.load.index
        self.line_set = net.line.index
        self.shunt_set = net.shunt.index
        self.trafo_set = net.trafo.index
        self.ext_grid_set = net.ext_grid.index
        self.slack_set = net.ext_grid.bus

        self.bus_gen_set = list(zip(self.generators.bus, self.gen_all_set))
        self.bus_demand_set = list(zip(net.load.bus, self.demand_set))
        self.bus_shunt_set = list(zip(net.shunt.bus, self.shunt_set))
        self.bus_ext_grid_set = list(zip(net.ext_grid.bus, self.ext_grid_set))

        self.bus_line_dict = dict(
            zip(list(zip(self.line_set, [1] * len(self.line_set))) + list(zip(self.line_set, [2] * len(self.line_set))),
                pd.concat([net.line.from_bus, net.line.to_bus])))
        self.bus_trafo_dict = dict(zip(list(zip(self.trafo_set, [1] * len(self.trafo_set))) + list(
            zip(self.trafo_set, [2] * len(self.trafo_set))), pd.concat([net.trafo.hv_bus, net.trafo.lv_bus])))

        # --- Param Data ---
        self.baseMVA = net.sn_mva

        self.PG_data = self.generators.p_mw * self.generators.scaling / self.baseMVA    # active power generation set points

        self.PD_data = self.net.load.p_mw * self.net.load.scaling / self.baseMVA

        self.GB_data = net.shunt.p_mw * net.shunt.step / self.baseMVA

        self.delta_eG_data = self.net.ext_grid.va_degree * pi/180

        # --- transformer ---
        self._calc_trafo_values()
        self.shift_rad_data = pd.Series(self.trafo_parameters["shift_rad"], self.trafo_set)
        self.tap_data = pd.Series(self.trafo_parameters["ratio"], self.trafo_set)

    def _calc_trafo_values(self):
        vn_trafo_hv, vn_trafo_lv, shift = self._calc_tap_shift()
        shift_rad = shift * pi/180
        ratio = self._calc_nominal_ratio_from_dataframe(vn_trafo_hv, vn_trafo_lv)
        r, x, y = self._calc_r_x_y_from_dataframe(vn_trafo_lv, self.net.bus.vn_kv[self.net.trafo.lv_bus].values)

        self.trafo_parameters = {'r': r, 'x': x, 'y': y, 'shift_rad': shift_rad, 'ratio': ratio}

    def _calc_tap_shift(self):
        """
        Adjust the nominal voltage vnh and vnl to the active tab position "tap_pos".
        If "side" is 1 (high-voltage side) the high voltage vnh is adjusted.
        If "side" is 2 (low-voltage side) the low voltage vnl is adjusted

        INPUT:
            **net** - The pandapower format network

            **trafo** (Dataframe) - The dataframe in pd_net["structure"]["trafo"]
            which contains transformer calculation values.

        OUTPUT:
            **vn_hv_kv** (1d array, float) - The adusted high voltages

            **vn_lv_kv** (1d array, float) - The adjusted low voltages

            **trafo_shift** (1d array, float) - phase shift angle

        """
        vnh = copy.deepcopy(self.net.trafo.vn_hv_kv.values)
        vnl = copy.deepcopy(self.net.trafo.vn_lv_kv.values)
        trafo_shift = self.net.trafo.shift_degree.values

        tap_pos = self.net.trafo.tap_pos
        tap_neutral = self.net.trafo.tap_neutral
        tap_diff = tap_pos - tap_neutral
        tap_phase_shifter = self.net.trafo.tap_phase_shifter
        tap_side = self.net.trafo.tap_side
        tap_step_percent = self.net.trafo.tap_step_percent
        tap_step_degree = self.net.trafo.tap_step_degree

        cos = lambda x: np.cos(np.deg2rad(x))
        sin = lambda x: np.sin(np.deg2rad(x))
        arctan = lambda x: np.rad2deg(np.arctan(x))

        for side, vn, direction in [("hv", vnh, 1), ("lv", vnl, -1)]:
            phase_shifters = tap_phase_shifter & (tap_side == side)
            tap_complex = np.isfinite(tap_step_percent) & np.isfinite(tap_pos) & (tap_side == side) & \
                          ~phase_shifters
            if tap_complex.any():
                tap_steps = tap_step_percent[tap_complex] * tap_diff[tap_complex] / 100
                tap_angles = (tap_step_degree[tap_complex]).fillna(0)
                u1 = vn[tap_complex]
                du = u1 * tap_steps.fillna(0)
                vn[tap_complex] = np.sqrt((u1 + du * cos(tap_angles)) ** 2 + (du * sin(tap_angles)) ** 2)
                trafo_shift[tap_complex] += (arctan(direction * du * sin(tap_angles) /
                                                    (u1 + du * cos(tap_angles))))
            if phase_shifters.any():
                degree_is_set = tap_step_degree[phase_shifters].fillna(0) != 0
                percent_is_set = tap_step_percent[phase_shifters].fillna(0) != 0
                if (degree_is_set & percent_is_set).any():
                    raise UserWarning(
                        "Both tap_step_degree and tap_step_percent set for ideal phase shifter")
                trafo_shift[phase_shifters] += np.where(
                    (degree_is_set),
                    (direction * tap_diff[phase_shifters] * tap_step_degree[phase_shifters]),
                    (direction * 2 * np.rad2deg(np.arcsin(tap_diff[phase_shifters] * \
                                                          tap_step_percent[phase_shifters] / 100 / 2)))
                )

        return vnh, vnl, trafo_shift

    def _calc_nominal_ratio_from_dataframe(self, vn_hv_kv, vn_lv_kv):
        """
        Calculates (Vectorized) the off nominal tap ratio::

                      (vn_hv_kv / vn_lv_kv) / (ub1_in_kv / ub2_in_kv)

        INPUT:
            **net** (Dataframe) - The net for which to calc the tap ratio.

            **vn_hv_kv** (1d array, float) - The adjusted nominal high voltages

            **vn_lv_kv** (1d array, float) - The adjusted nominal low voltages

        OUTPUT:
            **tab** (1d array, float) - The off-nominal tap ratio
        """
        # Calculating tab (trasformer off nominal turns ratio)
        tap_rat = vn_hv_kv / vn_lv_kv
        hv_bus = self.net.trafo.hv_bus
        lv_bus = self.net.trafo.lv_bus
        nom_rat = self.net.bus.vn_kv[hv_bus].values / self.net.bus.vn_kv[lv_bus].values
        return tap_rat / nom_rat

    def _calc_r_x_y_from_dataframe(self, vn_trafo_lv, vn_lv):
        trafo_model = self.net["_options"]["trafo_model"]

        r, x = self._calc_r_x_from_dataframe(vn_lv, vn_trafo_lv)

        y = self._calc_y_from_dataframe(vn_lv, vn_trafo_lv)

        if trafo_model == "pi":
            return r.values, x.values, y.values
        elif trafo_model == "t":
            return self._wye_delta(r.values, x.values, y.values)
        else:
            raise ValueError("Unkonwn Transformer Model %s - valid values ar 'pi' or 't'" % trafo_model)

    def _calc_r_x_from_dataframe(self, vn_lv, vn_trafo_lv):
        """
        Calculates (Vectorized) the resitance and reactance according to the
        transformer values
        """
        sn_mva = self.baseMVA
        parallel = self.net.trafo.parallel

        vk_percent = self.net.trafo.vk_percent
        vkr_percent = self.net.trafo.vkr_percent

        # adjust for low voltage side voltage converter:
        tap_lv = np.square(vn_trafo_lv / vn_lv) * sn_mva

        sn_trafo_mva = self.net.trafo.sn_mva
        z_sc = vk_percent / 100. / sn_trafo_mva * tap_lv
        r_sc = vkr_percent / 100. / sn_trafo_mva * tap_lv
        x_sc = np.sign(z_sc) * np.sqrt((z_sc ** 2 - r_sc ** 2).astype(float))
        return r_sc / parallel, x_sc / parallel

    def _calc_y_from_dataframe(self, vn_lv, vn_trafo_lv):
        """
        Calculate the subsceptance y from the transformer dataframe.

        INPUT:

            **trafo** (Dataframe) - The dataframe in net.trafo
            which contains transformer calculation values.

        OUTPUT:
            **subsceptance** (1d array, np.complex128) - The subsceptance in pu in

        """
        sn_mva = self.baseMVA

        baseR = np.square(vn_lv) / sn_mva
        vn_lv_kv = self.net.trafo.vn_lv_kv
        pfe = self.net.trafo.pfe_kw * 1e-3
        parallel = self.net.trafo.parallel

        ### Calculate subsceptance ###
        vnl_squared = vn_lv_kv ** 2
        b_real = pfe / vnl_squared * baseR
        i0 = self.net.trafo.i0_percent
        sn = self.net.trafo.sn_mva
        b_img = (i0 / 100. * sn) ** 2 - pfe ** 2

        b_img[b_img < 0] = 0
        b_img = np.sqrt(b_img) * baseR / vnl_squared
        y = b_real - 1j * b_img * np.sign(i0)
        return y / np.square(vn_trafo_lv / vn_lv_kv) * parallel

    def _wye_delta(self, r, x, y):
        """
        20.05.2016 added by Lothar Löwer

        Calculate transformer Pi-Data based on T-Data

        """
        y = y * -1j
        tidx = np.where(y != 0)
        za_star = (r[tidx] + x[tidx] * 1j) / 2
        zc_star = -1j / y[tidx]
        zSum_triangle = za_star * za_star + 2 * za_star * zc_star
        zab_triangle = zSum_triangle / zc_star
        zbc_triangle = zSum_triangle / za_star
        r[tidx] = zab_triangle.real
        x[tidx] = zab_triangle.imag
        y[tidx] = -2j / zbc_triangle
        y = y * 1j
        return r, x, y

    def create_model(self):
        self.model = ConcreteModel()

        # --- SETS ---
        self.model.B = Set(initialize=self.bus_set)
        self.model.eG = Set(initialize=self.ext_grid_set)  # external grids
        self.model.G = Set(initialize=self.gen_all_set)  # static generators and generators
        self.model.sG = Set(within=self.model.G, initialize=self.sgen_set)  # static generators
        self.model.gG = Set(within=self.model.G, initialize=self.gen_set)  # generators (not static)
        self.model.D = Set(initialize=self.demand_set)
        self.model.L = Set(initialize=self.line_set)
        self.model.SHUNT = Set(initialize=self.shunt_set)
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_set)

        # generators, buses, loads linked to each bus b
        self.model.Gbs = Set(within=self.model.B * self.model.G,
                             initialize=self.bus_gen_set)  # set of generator-bus mapping
        self.model.Dbs = Set(within=self.model.B * self.model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        self.model.SHUNTbs = Set(within=self.model.B * self.model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        self.model.eGbs = Set(within=self.model.B * self.model.eG,
                              initialize=self.bus_ext_grid_set)  # set of external grid-bus mapping

        # --- parameters ---
        # line matrix
        self.model.A = Param(self.model.L * self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF * self.model.LE,
                              initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # generation
        self.model.PG = Param(self.model.G, initialize=self.PG_data)

        # demand
        self.model.PD = Param(self.model.D, initialize=self.PD_data)

        self.model.VOLL = Param(self.model.D, within=Reals, initialize=10000)  # value of lost load

        # shunt
        self.model.GB = Param(self.model.SHUNT, within=Reals, initialize=self.GB_data)  # shunt conductance

        # trafo
        self.model.shift = Param(self.model.TRANSF, within=Reals,
                                 initialize=self.shift_rad_data)  # transformer phase shift in rad
        self.model.Tap = Param(self.model.TRANSF, within=Reals, initialize=self.tap_data)

        # external grid voltage angle
        self.model.delta_eG = Param(self.model.eG, within=Reals, initialize=self.delta_eG_data)

        # baseMVA of the net
        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)

        # --- variables ---
        self.model.delta = Var(self.model.B, domain=Reals, initialize=0.0,
                               bounds=(-pi, pi))  # voltage phase angle at bus b, rad
        self.model.pD = Var(self.model.D, domain=Reals)  # real power demand delivered
        self.model.pG = Var(self.model.G, domain=NonNegativeReals)  # real generator power
        self.model.peG = Var(self.model.eG, domain=Reals)  # real power injection from external grids
        self.model.pLfrom = Var(self.model.L, domain=Reals)  # real power injected at b onto line
        self.model.pLto = Var(self.model.L, domain=Reals)  # real power injected at b' onto line
        self.model.pThv = Var(self.model.TRANSF, domain=Reals)  # real power injected at b onto transformer
        self.model.pTlv = Var(self.model.TRANSF, domain=Reals)  # real power injected at b' onto transformer

        # --- generator power ---
        for g in self.model.G:
            self.model.pG[g].fix(self.model.PG[g])

        # --- demand ---
        for d in self.model.D:
            self.model.pD[d].fix(self.model.PD[d])

        # --- reference bus constraint ---
        for (b, g) in self.model.eGbs:
            self.model.delta[b].fix(self.model.delta_eG[g])
        # def ref_bus_def(model, b, g):
        #     return model.delta[b] == self.model.delta_eG[g]
        #
        # self.model.refbus_delta = Constraint(self.model.eGbs, rule=ref_bus_def)

    def solve(self, print_solver_output: bool = False, solver='ipopt', load_solutions: bool = True, mip_solver=None):
        optimizer = SolverFactory(solver)

        if solver == 'mindtpy':
            if mip_solver == 'gurobi':
                self.results = optimizer.solve(self.model, mip_solver='gurobi_persistent', nlp_solver='ipopt',
                                               tee=print_solver_output)
            else:
                self.results = optimizer.solve(self.model, mip_solver='glpk', nlp_solver='ipopt',
                                               tee=print_solver_output)
        else:
            self.results = optimizer.solve(self.model, load_solutions=load_solutions, tee=print_solver_output)
