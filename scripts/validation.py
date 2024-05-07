import pandas as pd
import simbench as sb
import pandapower as pp
import pickle
from tqdm import tqdm

from potpourri.models.AC import AC


if __name__ == "__main__":
    results_dir = '../potpourri/results'

    delta_keys = {'res_bus': ['vm_pu', 'va_degree'], 'res_line': ['pl_mw', 'ql_mvar']}

    nets = ['1-HV-mixed--0-sw', '1-HV-urban--0-sw', '1-MV-rural--0-sw', '1-MV-urban--0-sw', '1-LV-urban6--0-sw', '1-LV-rural1--0-sw']

    keys_flat = [item for sublist in delta_keys.values() for item in sublist]
    results_df = pd.DataFrame(columns=[keys_flat], index=nets)

    for net_name in tqdm(nets):
        print(net_name)
        net = sb.get_simbench_net(net_name)

        ac = AC(net)
        ac.solve()

        pp.runpp(net)
        pp.nets_equal(net, ac.net)

        pickle.dump(ac.net, open(results_dir + net_name + '_ac.pkl', 'wb'))

        for res_element, keys in delta_keys.items():
            for key in keys:
                delta = ac.net[res_element][key] - net[res_element][key]
                mean_delta = abs(delta).mean()
                results_df.loc[net_name, key] = mean_delta

    print(results_df.to_latex(float_format="{:.2e}".format,
                              header=['Amplitude [p.u.]', 'Winkel [°]', 'Wirkleistung [MW]', 'Blindleistung [MVar]']))
    # # Reset the index
    # results_df.reset_index(inplace=True)
    #
    # # Rename the index column to 'nets'
    # results_df.rename(columns={'index': 'nets'}, inplace=True)



