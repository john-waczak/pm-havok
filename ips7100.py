import numpy as np

pm_colors = [
    "#C3EE11",
    "#2CD387",
    "#2CCCD3",
    "#2C78D3",
    "#321CD3",
    "#852ED1",
    "#D12ECC",
]


pm_labels_units = {
    "0.1" : r"$\text{PM}_{0.1}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "0.3" : r"$\text{PM}_{0.3}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "0.5" : r"$\text{PM}_{0.5}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "1.0" : r"$\text{PM}_{1}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "2.5" : r"$\text{PM}_{2.5}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "5.0" : r"$\text{PM}_{5}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
    "10.0" : r"$\text{PM}_{10}$ $\left( \mu \text{g} \cdot \text{m}^{-3} \right)$",
}

pc_cols = ['pc0_1', 'pc0_3', 'pc0_5', 'pc1_0', 'pc2_5', 'pc5_0', 'pc10_0']
pm_cols = ['pm0_1', 'pm0_3', 'pm0_5', 'pm1_0', 'pm2_5', 'pm5_0', 'pm10_0']

bin_edges = [0.0, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0]



# assume df has pc columns
def get_conversion_factors(df):
        A = ((df.pm0_1[df.pm0_1.notna()] / df.pc0_1[df.pc0_1.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        B = (((df.pm0_3[df.pm0_3.notna()] - df.pm0_1[df.pm0_1.notna()] ) / df.pc0_3[df.pc0_3.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        C = (((df.pm0_5[df.pm0_5.notna()] - df.pm0_3[df.pm0_3.notna()] ) / df.pc0_5[df.pc0_5.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        D = (((df.pm1_0[df.pm1_0.notna()] - df.pm0_5[df.pm0_5.notna()] ) / df.pc1_0[df.pc1_0.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        E = (((df.pm2_5[df.pm2_5.notna()] - df.pm1_0[df.pm1_0.notna()] ) / df.pc2_5[df.pc2_5.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        F = (((df.pm5_0[df.pm5_0.notna()] - df.pm2_5[df.pm2_5.notna()] ) / df.pc5_0[df.pc5_0.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())
        G = (((df.pm10_0[df.pm10_0.notna()] - df.pm5_0[df.pm5_0.notna()] ) / df.pc10_0[df.pc10_0.notna()])
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .mean())

        return [A, B, C, D, E, F, G]


# def get_dry_bins(RH, kappa=0.62, eff=0.35, D_wet=np.array([0.0, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0])):
#     rh = 99.9 if RH == 100.0 else RH

#     if rh >= eff:
#         return D_wet / ((1 + kappa * rh / (100.0 - rh))**(1/3))
#     else:
#         return D_wet



def get_dry_bins(RH, kappa=0.62, eff=0.35):
    D_wet = np.tile(np.array(bin_edges), (len(RH), 1))
    correction = RH.copy()
    correction[RH == 100.0] = 99.99

    # create mask of values above efflorescence point
    eff_mask = correction >= eff

    # apply correction
    correction[eff_mask] = 1.0 / (1 + kappa * correction[eff_mask]/(100.0 - correction[eff_mask]))**(1/3)
    correction[~eff_mask] = 1.0
    
    D_wet = D_wet * correction[:, np.newaxis]

    return D_wet    



def pc_frac(bl, bh, Bl, Bh):
    return np.maximum((np.minimum(bh ,Bh) - np.maximum(bl, Bl))/(bh-bl), 0.0)

def pc_cor(pc_out, df, D_dry, Bl, Bh):    
    df[pc_out] = (
        df.pc0_1  * pc_frac(D_dry[:,0], D_dry[:,1], Bl, Bh)
        + df.pc0_3  * pc_frac(D_dry[:,1], D_dry[:,2], Bl, Bh)
        + df.pc0_5  * pc_frac(D_dry[:,2], D_dry[:,3], Bl, Bh)
        + df.pc1_0  * pc_frac(D_dry[:,3], D_dry[:,4], Bl, Bh)
        + df.pc2_5  * pc_frac(D_dry[:,4], D_dry[:,5], Bl, Bh)
        + df.pc5_0  * pc_frac(D_dry[:,5], D_dry[:,6], Bl, Bh)
        + df.pc10_0 * pc_frac(D_dry[:,6], D_dry[:,7], Bl, Bh)
    )



# assumes df has humidity, 
def correct_ips7100_for_humidity(df, kappa=0.62, eff=0.35, rh_col='humidity'):
    # get dry bis
    D_dry = get_dry_bins(df[rh_col].values, kappa=kappa, eff=eff)

    # now for each original pc, we need to 
    # redistribute counts based on 
    # the new bin edges    
    pc_cor('pc0_1_cor',  df, D_dry, bin_edges[0], bin_edges[1]) 
    pc_cor('pc0_3_cor',  df, D_dry, bin_edges[1], bin_edges[2]) 
    pc_cor('pc0_5_cor',  df, D_dry, bin_edges[2], bin_edges[3]) 
    pc_cor('pc1_0_cor',  df, D_dry, bin_edges[3], bin_edges[4]) 
    pc_cor('pc2_5_cor',  df, D_dry, bin_edges[4], bin_edges[5]) 
    pc_cor('pc5_0_cor',  df, D_dry, bin_edges[5], bin_edges[6]) 
    pc_cor('pc10_0_cor', df, D_dry, bin_edges[6], bin_edges[7]) 

    # compute conversion constants    
    coeffs = get_conversion_factors(df)

    # get correct PM values
    df['pm0_1_cor']  = coeffs[0] * df.pc0_1_cor
    df['pm0_3_cor']  = coeffs[1] * df.pc0_3_cor + df.pm0_1_cor
    df['pm0_5_cor']  = coeffs[2] * df.pc0_5_cor + df.pm0_3_cor
    df['pm1_0_cor']  = coeffs[3] * df.pc1_0_cor + df.pm0_5_cor
    df['pm2_5_cor']  = coeffs[4] * df.pc2_5_cor + df.pm1_0_cor
    df['pm5_0_cor']  = coeffs[5] * df.pc5_0_cor + df.pm2_5_cor
    df['pm10_0_cor'] = coeffs[6] * df.pc10_0_cor + df.pm5_0_cor

