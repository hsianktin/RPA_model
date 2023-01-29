# estimate the chemical constants based on fitted parameters
using DataFrames,CSV
using Pkg
Pkg.add("Unitful")
using Unitful


para_df = CSV.read("figs/para_fitted.csv",DataFrame)

cₒ = 1e-11u"M" # M

k₊ = para_df.k_on ./ cₒ /1u"s" # M⁻¹s⁻¹
k₋ = para_df.k_off /1u"s"  # s⁻¹
Kd₂₀ = k₋ ./ k₊ # M
Kd₃₀ = Kd₂₀ .* para_df.v_open ./ (para_df.v_open .+ para_df.v_close ) # M

kinetics_df = DataFrame(
    salt_concentration = para_df.exp_label,
    k₊₂₀ = k₊,
    k₋₂₀ = k₋,
    Kd₂₀ = Kd₂₀,
    Kd₃₀ = Kd₃₀,
)

CSV.write("figs/sources/kinetics.csv",kinetics_df)