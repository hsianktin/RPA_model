using Base: Float64
using DataFrames
using CSV
datapath = "./data_simu"
dataname(id) = "$datapath/diffusion_report_$(id).csv"
ids = [i for i in 1:18]
summarizing_table = DataFrame(wt_15_1_20=Float64[],wt_15_1_30=Float64[],wt_15_2_20=Float64[],wt_15_2_30=Float64[],wt_150_1_20=Float64[],wt_150_1_30=Float64[],wt_150_2_20=Float64[],wt_150_2_30=Float64[],dwh_150_1_20=Float64[],dwh_150_1_30=Float64[],dwh_150_2_20=Float64[],dwh_150_2_30=Float64[],df_150_1_20=Float64[],df_150_1_30=Float64[],df_150_2_20=Float64[],df_150_2_30=Float64[],wt_15_1_q=Float64[],wt_15_2_q=Float64[],wt_150_1_q=Float64[],wt_150_2_q=Float64[],dwh_150_1_q=Float64[],dwh_150_2_q=Float64[],df_150_1_q=Float64[],df_150_2_q=Float64[])
for id in ids
    df = CSV.read(dataname(id),DataFrame)
    data = Array{Float64,1}()
    list = [1,2,3,4,6,7,8,9]
    for line in list
        push!(data,df[line,4])
        push!(data,df[line,5])
    end
    for line in list
        push!(data,df[line,3])
    end
    push!(summarizing_table,data)
end
CSV.write("./data_simu/summarizing_fret_diffusion.csv",summarizing_table)