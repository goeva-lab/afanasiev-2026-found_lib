# %%
import altair as alt
import anndata as ad
from altair.theme import register
from found import HiDDENg, pl
from scipy import sparse as sp


def set_wh(h: int, w: int):
    register("cust", enable=True)(lambda: {"height": h, "width": w})


alt.renderers.enable("svg")

# %%
adata = ad.concat([ad.read_h5ad(f"./.cache/vignette/{mol}.h5ad") for mol in ["pbs", "il15"]])
adata.X = sp.csr_array(adata.X)

# %%
adata.obs["cytokine"] = adata.obs["cytokine"].cat.reorder_categories(["PBS", "IL-15"], ordered=False)  # ty:ignore[unresolved-attribute]

# %%
p_hat, labs = HiDDENg(adata, k=50, group_by="cell_type", cond_col="cytokine", control_val="PBS", X=adata.X)
adata.obs["HiDDEN_p_hat"] = p_hat
adata.obs["HiDDEN_labs"] = labs


# %%
plt = pl.PlotHiDDENOutput(adata, p_hat, labs)
# %%
set_wh(200, 200)
# %%
plt.reg_vln("cytokine")
# %%
set_wh(200, 245)
# %%
plt.bin_bar("cytokine", "PBS", "cell_type", vertical=False)
