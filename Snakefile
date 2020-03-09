import itertools

nsim = [100]
ngene = [10000]
pdeg = [0.05, 0.25]
assign = [[1/2, 1/2], [0.7, 0.3], [0.9, 0.1], [1.0, 0.0]]
foldchange = [[4, 4]]
replicates = [[3, 3]]
names = ["nsim", 'ngene', 'pdeg', 'assign', 'foldchange', 'replicates']

sim_params = {f"condition{i+1:02}":
 {k:l for k, l in zip(names, x)}
 for i, x in enumerate(itertools.product(nsim, ngene, pdeg, assign, foldchange, replicates))}

sim_params2 = {f"condition{i+1:02}":
 {k:l for k, l in zip(names, x)}
 for i, x in enumerate(itertools.product(nsim, ngene, pdeg, assign[0:-1], foldchange, replicates))}

assign = [[1/3, 1/3, 1/3], [0.6, 0.2, 0.2], [0.8, 0.1, 0.1], [0.5, 0.5, 0.0], [1.0, 0.0, 0.0]]
foldchange = [[4, 4, 4]]
replicates = [[3, 3, 3]]

sim_params3 = {f"condition{i+1:02}":
 {k:l for k, l in zip(names, x)}
 for i, x in enumerate(itertools.product(nsim, ngene, pdeg, assign, foldchange, replicates))}

rule all:
    input: 
        #expand("results/two_groups/{condition}/mean_fdr_q.rds", condition = sim_params.keys()),
        #expand("results/two_groups_tcc/{condition}/mean_fdr_q.rds", condition = sim_params2.keys()),
        #expand("results/three_groups/{condition}/mean_fdr_q.rds", condition = sim_params3.keys()),
        #expand("results/fc_two_groups/{condition}/mean_fdr_q.rds", condition = sim_params.keys()),
        #expand("results/two_groups/{condition}/summarize_result.rds", condition = sim_params.keys()),
        #expand("results/two_groups_tcc/{condition}/summarize_result.rds", condition = sim_params2.keys()),
        #expand("results/three_groups/{condition}/summarize_result.rds", condition = sim_params3.keys()),
        #expand("results/fc_two_groups/{condition}/summarize_result.rds", condition = sim_params.keys()),
        #expand("results/fc_two_groups_nstart/{condition}/summarize_result.rds", condition = sim_params.keys()),
        #expand("results/three_groups/{condition}/summarize_result.rds", condition = sim_params3.keys()),
        #expand("results/fc_two_groups_tcc/{condition}/summarize_result.rds", condition = sim_params.keys()),
        #expand("results/fc_two_groups_nstart2/{condition}/summarize_result.rds", condition = sim_params.keys()),
        #"results/all_metrics.rds",
        #"results/fdr_q.html",
        #"results/real_data/result.rds",
        #"results/two_groups_simseq/all_results.rds",
        "results/three_groups_simseq/all_results.rds",
        
rule bench_2groups:
    output: 
        "results/two_groups/{condition}/all_results.rds",
        "results/two_groups/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params[wildcards.condition]
    threads:
        12
    log:
        "logs/two_groups/{condition}/stdout.log",
        "logs/two_groups/{condition}/stderr.log",
    script:
        "scripts/bench_2groups.R"
        
rule bench_2groups_simseq:
    output: 
        "results/two_groups_simseq/all_results.rds",
        #"results/two_groups_simseq/params.json"
    threads:
        12
    log:
        "logs/two_groups_simseq/stdout.log",
        "logs/two_groups_simseq/stderr.log",
    script:
        "scripts/simseq_2groups.R"

rule bench_3groups_simseq:
    output: 
        "results/three_groups_simseq/all_results.rds",
        #"results/three_groups_simseq/params.json"
    threads:
        12
    log:
        "logs/three_groups_simseq/stdout.log",
        "logs/three_groups_simseq/stderr.log",
    script:
        "scripts/simseq_3groups.R"

rule bench_2groups_mbcluster_with_tcc:
    output: 
        "results/two_groups_tcc/{condition}/all_results.rds",
        "results/two_groups_tcc/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params2[wildcards.condition]
    threads:
        12
    log:
        "logs/two_groups_tcc/{condition}/stdout.log",
        "logs/two_groups_tcc/{condition}/stderr.log",
    script:
        "scripts/tcc_mbcluster.R"

#rule bench_2groups_mbcluster_with_zero_size:
#    output: 
#        "results/two_groups_zero/{condition}/all_results.rds",
#        "results/two_groups_zero/{condition}/params.json"
#    params:
#        lambda wildcards, output: sim_params2[wildcards.condition]
#    threads:
#        12
#    log:
#        "logs/two_groups_tcc/{condition}/stdout.log",
#        "logs/two_groups_tcc/{condition}/stderr.log",
#    script:
#        "scripts/tcc_mbcluster.R"

rule bench_3groups:
    output: 
        "results/three_groups/{condition}/all_results.rds",
        "results/three_groups/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params3[wildcards.condition]
    threads:
        12
    log:
        "logs/three_groups/{condition}/stdout.log",
        "logs/three_groups/{condition}/stderr.log",
    script:
        "scripts/bench_3groups.R"

rule bench_fc_2groups:
    output: 
        "results/fc_two_groups/{condition}/all_results.rds",
        "results/fc_two_groups/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params[wildcards.condition]
    threads:
        12
    log:
        "logs/fc_two_groups/{condition}/stdout.log",
        "logs/fc_two_groups/{condition}/stderr.log",
    script:
        "scripts/bench_fc_2groups.R"

rule bench_fc_2groups_nstart:
    output: 
        "results/fc_two_groups_nstart/{condition}/all_results.rds",
        "results/fc_two_groups_nstart/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params[wildcards.condition]
    threads:
        12
    log:
        "logs/fc_two_groups_nstart/{condition}/stdout.log",
        "logs/fc_two_groups_nstart/{condition}/stderr.log",
    script:
        "scripts/bench_fc_2groups_nstart.R"

rule bench_fc_2groups_nstart2:
    output: 
        "results/fc_two_groups_nstart2/{condition}/all_results.rds",
        "results/fc_two_groups_nstart2/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params[wildcards.condition]
    threads:
        12
    log:
        "logs/fc_two_groups_nstart2/{condition}/stdout.log",
        "logs/fc_two_groups_nstart2/{condition}/stderr.log",
    script:
        "scripts/bench_fc_2groups_nstart2.R"

rule bench_fc_2groups_tcc:
    output: 
        "results/fc_two_groups_tcc/{condition}/all_results.rds",
        "results/fc_two_groups_tcc/{condition}/params.json"
    params:
        lambda wildcards, output: sim_params[wildcards.condition]
    threads:
        12
    log:
        "logs/fc_two_groups_tcc/{condition}/stdout.log",
        "logs/fc_two_groups_tcc/{condition}/stderr.log",
    script:
        "scripts/bench_fc_2groups_tcc.R"

rule summarize_results:
    input: 
        "results/{design}/{condition}/all_results.rds",
        "results/{design}/{condition}/params.json",
    output:
        "results/{design}/{condition}/summarize_result.rds",
    threads: 3
    log:
        "logs/summarize/{design}/{condition}/stdout.log",
        "logs/summarize/{design}/{condition}/stderr.log",
    script:
        "scripts/summarize_results.R"

rule get_mean_fdr_q:
    input: 
        "results/{design}/{condition}/summarize_result.rds",
    output:
        "results/{design}/{condition}/mean_fdr_q.rds"
    
    script:
        "scripts/get_mean_fdr_q.R"

#rule calc_metrics:
#    input:
#        expand("results/two_groups/{condition}/summarize_result.rds", condition = sim_params.keys()),
#        expand("results/two_groups_tcc/{condition}/summarize_result.rds", condition = sim_params2.keys()),
#        expand("results/three_groups/{condition}/summarize_result.rds", condition = sim_params3.keys()),
#        expand("results/fc_two_groups/{condition}/summarize_result.rds", condition = sim_params.keys()),
#    output:
#        "results/all_metrics.rds",
#    script:
#        "scripts/calc_metrics.R"

rule report_fdr_q:
    input: 
        expand("results/two_groups/{condition}/mean_fdr_q.rds", condition = sim_params.keys()),
        #expand("results/two_groups_tcc/{condition}/mean_fdr_q.rds", condition = sim_params2.keys()),
        #expand("results/three_groups/{condition}/mean_fdr_q.rds", condition = sim_params3.keys()),
        #expand("results/fc_two_groups/{condition}/mean_fdr_q.rds", condition = sim_params.keys()),
    output:
        "results/fdr_q.html"
    
    script:
        "scripts/report_fdr_q.Rmd"

rule real_data:
    output:
        "results/real_data/result.rds"
    threads:12
    script:
        "scripts/real_data.R"
