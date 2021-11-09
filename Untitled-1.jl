exp_label = "wt_15mM_salt"
simu_labels = ["further", "specific_nd"]
for simu_label in simu_labels
try
open("$simupath/rsa_state_transition_$(exp_label)_$(simu_label).csv","r") do input
        open("$simupath/rsa_state_transition_$(exp_label).csv","a") do output
            n = countlines(input)
            # record = 0
            seekstart(input)
            for k in 1:n
                line = readline(input)
                println(output,line)
                # record = record + 1
            end
            # println(record)
        end
    end
    catch
    end
end