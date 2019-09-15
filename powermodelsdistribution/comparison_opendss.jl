import PowerModelsDistribution
PMD = PowerModelsDistribution
import DistributionTestCases
DTC = DistributionTestCases

import JuMP
import Ipopt
ipopt_solver = JuMP.with_optimizer(Ipopt.Optimizer, tol=1E-10)

# name of cases
cases = ["IEEE13", "IEEE34", "IEEE123", "LVTestCase t=500", "LVTestCase t=1000"]
# path for each case
paths = Dict(
        "IEEE13"            => DTC.CASE_PATH["IEEE13"],
        "IEEE34"            => DTC.CASE_PATH["IEEE34"],
        "IEEE123"           => DTC.CASE_PATH["IEEE123"],
        "LVTestCase t=500"  => DTC.CASE_PATH["LVTestCase"][500],
        "LVTestCase t=1000" => DTC.CASE_PATH["LVTestCase"][1000],
)

# indicate for each case which buses should be compared line-to-line instead of line-to-neutral
buses_compare_ll = Dict([(case, []) for case in cases])
buses_compare_ll["IEEE123"] = ["610"]

# gather the results
results = Dict()
for case in cases
        path = paths[case]
        results[case] = Dict()
        sol_dss = DTC.get_soldss_opendssdirect(path, tolerance=1E-6)
        data_pmd = PMD.parse_file(path)

        sol_pmd = PMD.run_ac_mc_pf_lm(data_pmd, ipopt_solver)

        results[case]["δ"] = DTC.compare_dss_to_pmd(sol_dss, sol_pmd, data_pmd, buses_compare_ll=buses_compare_ll[case])

        vm_all = [bus["vm"][c] for (_, bus) in sol_pmd["solution"]["bus"] for c in 1:3]
        results[case]["vm_max"] = maximum(vm_all)
        results[case]["vm_min"] = minimum(vm_all)
end
results
