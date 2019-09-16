import PowerModels
PMs = PowerModels
import PowerModelsDistribution
PMD = PowerModelsDistribution
import DistributionTestCases
DTC = DistributionTestCases

import JuMP
import Ipopt
ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer, tol=1E-8)
import Mosek, MosekTools
mosek_solver = JuMP.with_optimizer(Mosek.Optimizer)

# name of cases
#cases = ["IEEE13", "IEEE34", "IEEE123", "LVTestCase t=500", "LVTestCase t=1000"]
cases = ["IEEE13", "IEEE34", "LVTestCase t=500", "LVTestCase t=1000"]

# path for each case
paths = Dict(
        "IEEE13"            => DTC.CASE_PATH["IEEE13"],
        "IEEE34"            => DTC.CASE_PATH["IEEE34"],
        "IEEE123"           => DTC.CASE_PATH["IEEE123"],
        "LVTestCase t=500"  => DTC.CASE_PATH["LVTestCase"][500],
        "LVTestCase t=1000" => DTC.CASE_PATH["LVTestCase"][1000],
)

# scale down load to keep vm in range [0.9, 1.1]
load_scale = Dict([(case, 0.5) for case in cases])

function get_bound_s(data_pmd)
        pd_total = sum(vcat([load["pd"] for (_, load) in data_pmd["load"]]...))
        qd_total = sum(vcat([load["qd"] for (_, load) in data_pmd["load"]]...))
        sd_total = sqrt(pd_total^2+qd_total^2)
        return sd_total*2
end

# function to add OPF specific data to the data model
function add_opf_data!(data_pmd; vmin=0.9, vmax=1.1, sbound=get_bound_s(data_pmd))

        for (_, br) in data_pmd["branch"]
                br["rate_a"] = PMs.MultiConductorVector(ones(3)*sbound)
        end

        for (_, bus) in data_pmd["bus"]
                if bus["name"]=="virtual_sourcebus"
                        vm = bus["vm"].values
                        bus["vmin"] = PMs.MultiConductorVector(vm*1)
                        bus["vmax"] = PMs.MultiConductorVector(vm*1)
                else
                        bus["vmin"] = PMs.MultiConductorVector(ones(3)*vmin)
                        bus["vmax"] = PMs.MultiConductorVector(ones(3)*vmax)
                end
        end

        sourcegen = [gen for (_, gen) in data_pmd["gen"] if gen["name"]=="sourcegen"][1]
        sourcegen["pmax"] =  PMs.MultiConductorVector(ones(3)*sbound)
        sourcegen["pmin"] = -PMs.MultiConductorVector(ones(3)*sbound)
        sourcegen["qmax"] =  PMs.MultiConductorVector(ones(3)*sbound)
        sourcegen["qmin"] = -PMs.MultiConductorVector(ones(3)*sbound)
end
function get_obj(pm)
        try
                return JuMP.objective_value(pm.model)
        catch
                obj = 0
                for id in PMs.ids(pm, :gen)
                        pg = JuMP.value.([PMs.var(pm, pm.cnw, c, :pg, id) for c in 1:3])
                        cost = PMs.ref(pm, :gen, id, "cost")
                        obj += sum(cost.*pg)
                end
                return obj
        end
end

# gather the results
results = Dict()
for case in cases
        results[case] = Dict()

        path = paths[case]
        data_pmd = PMD.parse_file(path)
        DTC.simplify_feeder!(data_pmd, load_scale=load_scale[case])
        add_opf_data!(data_pmd)

        # exact formulations
        exact_form = [PMD.ACPPowerModel, PMD.ACRPowerModel]
        exact_label = ["AC-POLAR-NLP", "AC-RECT-NLP"]

        for (i, form) in enumerate(exact_form)
                pm = PMs.build_model(data_pmd, form, PMD.post_mc_opf_lm; ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
                sol_pmd = PMs.optimize_model!(pm, ipopt_solver)

                results_exact = Dict("objective"=>NaN, "termination_status"=>sol_pmd["termination_status"])
                results[case][exact_label[i]] = results_exact
                if sol_pmd["termination_status"]==PMD.LOCALLY_SOLVED
                        obj_relax = get_obj(pm)
                        results_exact["objective"] = get_obj(pm)
                        results_exact["solve_time"] = sol_pmd["solve_time"]
                end
        end

        # relaxations
        relax_form = [PMD.SDPUBFPowerModel, PMD.SOCConicUBFPowerModel]
        relax_label = ["BFM-SDP", "BFM-SOC"]
        obj_exact = results[case]["AC-POLAR-NLP"]["objective"]

        for (i, form) in enumerate(relax_form)
                pm = PMs.build_model(data_pmd, form, PMD.post_mc_opf_bf; ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
                sol_pmd = PMs.optimize_model!(pm, mosek_solver)
                results_relax = Dict("objective"=>NaN, "termination_status"=>sol_pmd["termination_status"])
                results[case][relax_label[i]] = results_relax
                if sol_pmd["termination_status"]==PMD.OPTIMAL
                        obj_relax = get_obj(pm)
                        results_relax["objective"] = obj_relax
                        results_relax["solve_time"] = sol_pmd["solve_time"]
                        results_relax["gap"] = (obj_exact-obj_relax)/obj_exact
                end
        end
end
results
