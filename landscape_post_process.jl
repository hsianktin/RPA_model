using DataFrames, CSV
# parse L from argument input 
L = parse(Int, ARGS[1])


df = CSV.read("data_simu/rsa_plot__landscape_L_$(L).csv.gz", DataFrame)
cols = [
    "k_on",
    "k_off",
    "v_open",
    "v_close",
    "fold",
    "L",
    "T1",
    "T2",
    "D1",
    "Id",
]
dt = 5
ts = [string(i) for i in 1800:dt:2400] 
selected_cols = [cols; ts]

df = df[:, selected_cols]

CSV.write("data_simu/rsa_plot_landscape_selected_L_$(L).csv.gz", df; compress=true)