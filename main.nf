// Perform standard QC steps
process MISSINGQC {
    
    input:
    val sample

    output:
    tuple path("*.smiss"), path("*.vmiss")
    tuple path("${sample}_postQC.bed"), path("${sample}_postQC.bim"), path("${sample}_postQC.fam")

    script:
    """
    plink2 --bfile ${params.dataset}/${sample} --missing --out ${sample}
    plink2 --bfile ${params.dataset}/${sample} --geno ${task.ext.relaxed_thres} --make-bed --out ${sample}_filt1
    plink2 --bfile ${sample}_filt1 --mind ${task.ext.relaxed_thres} --make-bed --out ${sample}_filt2
    plink2 --bfile ${sample}_filt2 --geno ${task.ext.stringent_thres} --make-bed --out ${sample}_filt3
    plink2 --bfile ${sample}_filt3 --mind ${task.ext.stringent_thres} --make-bed --out ${sample}_postQC
    """
}

// Plot missingness
process PLOTMISSING {

    input:
    tuple path(smiss), path(vmiss)

    output:
    path("$task.ext.outfile")

    script:
    """
    python ${params.script}/plot_missingness.py $smiss $vmiss $task.ext.outfile
    """
}

// Create ethnic group tables
process CREATE_ETHNIC_TABLES {
    input:
    val ethnic_file

    output:
    path "*.txt"

    script:
    """
    python ${params.script}/split_ethnicity.py ${ethnic_file}
    """
}

// Split samples based on ethnicity
process SPLITSAMPLES {
    input:
    tuple val(base_path), path(ethnic_file)

    output:
    tuple path("*_${ethnic_file.baseName}.bed"), 
          path("*_${ethnic_file.baseName}.bim"), 
          path("*_${ethnic_file.baseName}.fam")

    script:
    def new_base_path = base_path.tokenize('/').last()
    """
    plink2 --bfile $base_path --keep-fam $ethnic_file --make-bed --out ${new_base_path}_${ethnic_file.baseName}
    """
}

// Conduct GWAS analysis for each ethnic group
process GWAS {
    input:
    tuple path(bedfile), path(bimfile), path(famfile)

    output:
    path "*.linear"

    script:
    def base_path = bedfile.toString().replaceAll(/\.bed$/, '')
    """
    plink2 --bfile $base_path --glm allow-no-covars --pheno ${params.pheno} --out $base_path
    """
}

// Identify common variants and output the associated QQplot and Manhattan plot
process COMMONVARIANT {
    
    input:
    tuple val(var_string), path(files)

    output:
    path "${var_string}.csv" optional true
    path "*.png" optional true

    script:
    """
    python ${params.script}/common_variant.py -i $files -p $task.ext.pval -o ${var_string}.csv -w $task.ext.weird
    """
}


// Run workflow
workflow {
    // Run QC and plot missingness
    (missing_ch, filtered_ch) = MISSINGQC(params.sample)
    PLOTMISSING(missing_ch)
    
    // Get the base path for the genotype files in the tuple
    base_path_ch = filtered_ch.map { tuple ->
       tuple[1].toString()replaceAll(/\.(bed|bim|fam)$/, '')
       }
    
    // Create ethnic group tables and combine with base path
    ethnic_ch = CREATE_ETHNIC_TABLES(params.ethnic)
        .flatten()
    combined_ch = base_path_ch
        .combine(ethnic_ch)
        .map { base_path, ethnic -> [base_path, ethnic] }

    // Split genotype samples based on ethnicity and perform GWAS
    ethnic_snps_ch = SPLITSAMPLES(combined_ch)
    gwas_ch = GWAS(ethnic_snps_ch)

    // Identify common variants from GWAS results
    gwas_by_ethnic_ch = gwas_ch
        .flatten()
        .map { file ->
            tuple(file.name.replaceAll(/_eth_[^.]+\./, '.'), file)
        }
        .groupTuple()
    (csv_ch, png_ch) = COMMONVARIANT(gwas_by_ethnic_ch)
}