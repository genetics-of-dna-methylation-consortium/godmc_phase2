#!/usr/bin/env bash
source "${scripts_directory:-.}/resources/genetics/ld_workflow.sh"
source "${scripts_directory:-.}/resources/genetics/ld_pack.sh"

check_results_01 () {

	if [ -f "${cohort_descriptives}" ]; then
		echo "Cohort descriptives file present"
	else
		echo "Cohort descriptives file absent. Please re-run."
		exit 1
	fi

	#if [ -f "${methylation_summary}" ]; then
	#	echo "Methylation summary file present"
	#else
	#	echo "Methylation summary file absent. Please re-run."
	#	exit 1
	#fi

}


check_results_02 () {

	if [ -f "${allele_ref}" ]; then
		echo "Allele reference file present"
	else
		echo "Problem: Allele reference file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/pcaplot.pdf" ]; then
		echo "PCA plot present"
	else
		echo "Problem: PCA plot is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/easyQC_hrc.multi.AFCHECK.png" ]; then
		echo "easyQC plot present"
	else
		echo "Problem: easyQC plot is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.afreq.gz" ]; then
		echo "Allele frequency file present"
	else
		echo "Problem: Allele frequency file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.hardy.gz" ]; then
		echo "HWE file present"
	else
		echo "Problem: HWE file is absent"
		exit 1
	fi

	if [ -f "${section_02_dir}/data.info.gz" ]; then
		echo "Imputation quality file present"
	else
		echo "Problem: Imputation quality file is absent"
		exit 1
	fi

}

check_results_03a () {

	if [ -f "${section_03_dir}/methylation_summary.RData" ]; then
		echo "Methylation_summary.RData is present"
	else
		echo "Problem: methylation_summary.RData is absent"
		exit 1
	fi

    if [ -f "${section_03_dir}/cohort_descriptives_commonids.RData" ]; then
		echo "cohort_descriptives_commonids.RData is present"
	else
		echo "Problem: cohort_descriptives_commonids.RData is absent"
		exit 1
	fi

	
	if [ -f "${section_03_dir}/cellcounts_summary.txt" ]; then
		echo "Summary statistics of cell counts are present"
	else
		echo "Problem: summary statistics of cell counts are absent"
		exit 1
	fi

	if [ -f "${section_03_dir}/cellcounts_plot.pdf" ]; then
		echo "Plots of cell counts are present"
	else
		echo "Problem: plots of cell counts are absent"
		exit 1
	fi


	if [ "${measured_cellcounts}" != "NULL" ];then
		if [  -f "${section_03_dir}/cor_plot.pdf" ]; then
			echo "Correlation plot of observed vs predicted cell counts is present"
		else
			echo "Problem: correlation plot of observed vs predicted cell counts is absent"
			exit 1
        fi        

		if [  -f "${section_03_dir}/cor_matrix.txt" ]; then
			echo "Correlation matrix of observed vs predicted cell counts is present"
		else
			echo "Problem: correlation matrix of observed vs predicted cell counts is absent"
			exit 1
		fi
	else
		echo "Message: since measured_cellcounts are not provided, there is no output for cor_plot.pdf and cor_matrix.txt for observed vs predicted cell counts."

	fi

	if [ -f "${smoking_pred_plot}" ]; then
		echo "Smoking prediction plot is present"
	else
		echo "Problem: Smoking prediction plot file not present"
		exit 1
	fi

	if [ -f "${section_03_dir}/age_prediction.pdf" ]; then
		echo "The correlation plot between predicted and actual ages is present"
	else
		echo "Problem: The correlation plot between predicted and actual ages is absent"
		exit 1
	fi

	same_age_check=$(awk 'NR>1 {print $3}' ${covariates} | sort -n | uniq | wc -l)
	if [ "${same_age_check}" -eq 1 ]; then
		echo "All individuals have the same age. Skipping age prediction correlation matrix and statistics."
	else
		if [ -f "${section_03_dir}/age_prediction_correlation.png" ]; then
			echo "The matrix correlation plot among predicted age, age acceleration residuals, and chronological age is present"
		else
			echo "Problem: The matrix correlation plot of predicted age is absent"
			exit 1
		fi
	fi

	if [ -f "${section_03_dir}/age_prediction_stats.csv" ]; then
		echo "The statistical table for each clock and their age acceleration modules is present"
	else
		echo "Problem: The statistical table for each clock and their age acceleration modules is absent"
		exit 1
	fi

	if [ -f "${section_03_dir}/age_prediction_stats_corrsd.csv" ]; then
		echo "The statistical table for each pair of comparisons in the matrix plot of aging is present"
	else
		echo "Problem: The statistical table for each pair of comparisons in the matrix plot of aging is absent"
		exit 1
	fi

}

check_results_03d () {
    	echo "The number of methylation files may varied across cohorts. For more details, please check the Wiki. Please ensure you have seen the scripts 03a-03d run successfully from log files."
}

check_results_03 () {

	check_results_03a

 	check_results_03d

 	if [ -f "${section_03_dir}/positive_control_transformed_${positive_control_cpg}.PHENO1.glm.linear.gz" ]; then
		echo "transformed mQTL analysis positive control results present"
	else
		echo "transformed mQTL analysis positive control results file not present"
		exit 1
	fi

	if [ -f "${section_03_dir}/positive_control_transformed_${positive_control_cpg}_manhattan.pdf" ]; then
		echo "transformed mQTL analysis positive control Manhattan plot present"
	else
		echo "transformed mQTL analysis positive control Manhattan plot file not present"
		exit 1
	fi

		if [ -f "${section_03_dir}/positive_control_transformed_${positive_control_cpg}_nocisChr_manhattan.pdf" ]; then
		echo "transformed mQTL analysis positive control Manhattan plot present"
	else
		echo "transformed mQTL analysis positive control Manhattan plot file not present"
		exit 1
	fi

	if [ -f "${section_03_dir}/positive_control_transformed_${positive_control_cpg}_qqplot.jpeg" ]; then
		echo "transformed mQTL analysis positive control QQ plot present"
	else
		echo "transformed mQTL analysis positive control QQ plot file not present"
		exit 1
	fi

		if [ -f "${section_03_dir}/positive_control_transformed_${positive_control_cpg}_nocisChr_qqplot.jpeg" ]; then
		echo "transformed mQTL analysis positive control QQ plot present"
	else
		echo "transformed mQTL analysis positive control QQ plot file not present"
		exit 1
	fi

	if [ -f "${section_03_dir}/positive_control_untransformed_${positive_control_cpg}.PHENO1.glm.linear.gz" ]; then
		echo "untransformed analysis positive control results present"
	else
		echo "untransformed analysis positive control results file not present"
		exit 1
	fi
	
	if [ -f "${section_03_dir}/positive_control_untransformed_${positive_control_cpg}_manhattan.pdf" ]; then
		echo "untransformed analysis positive control Manhattan plot present"
	else
		echo "untransformed analysis positive control Manhattan plot file not present"
		exit 1
	fi

		if [ -f "${section_03_dir}/positive_control_untransformed_${positive_control_cpg}_nocisChr_manhattan.pdf" ]; then
		echo "untransformed analysis positive control Manhattan plot present"
	else
		echo "untransformed analysis positive control Manhattan plot file not present"
		exit 1
	fi
	
	if [ -f "${section_03_dir}/positive_control_untransformed_${positive_control_cpg}_qqplot.jpeg" ]; then
		echo "untransformed analysis positive control QQ plot present"
	else
		echo "untransformed analysis positive control QQ plot file not present"
		exit 1
	fi

		if [ -f "${section_03_dir}/positive_control_untransformed_${positive_control_cpg}_nocisChr_qqplot.jpeg" ]; then
		echo "untransformed analysis positive control QQ plot present"
	else
		echo "untransformed analysis positive control QQ plot file not present"
		exit 1
	fi
	
}

check_results_04 () {

	if [ -f "${home_directory}/results/04/meta_inputs/part_dev/${study_name}_a_cov.npy" ]; then
		echo "${study_name}_a_cov.npy present"
	else
		echo "${study_name}_a_cov.npy absent. Please re-run."
		exit 1
	fi

	if [ -f "${home_directory}/results/04/meta_inputs/part_dev/${study_name}_C.npy" ]; then
		echo "${study_name}_C.npy present"
	else
		echo "${study_name}_C.npy absent. Please re-run."
		exit 1
	fi

	if [ -f "${home_directory}/results/04/meta_inputs/part_dev/${study_name}_b_cov.npy" ]; then
		echo "${study_name}_b_cov.npy present"
	else
		echo "${study_name}_b_cov.npy absent. Please re-run."
		exit 1
	fi

	if [ -f "${home_directory}/results/04/meta_inputs/part_dev/${study_name}_a_test.npy" ]; then
		echo "${study_name}_a_test.npy present"
	else
		echo "${study_name}_a_test.npy absent. Please re-run."
		exit 1
	fi

	if [ -f "${home_directory}/results/04/meta_inputs/part_dev/${study_name}_metadata.npy" ]; then
		echo "${study_name}_metadata.npy present"
	else
		echo "${study_name}_metadata.npy absent. Please re-run."
		exit 1
	fi

    if [ -f "${home_directory}/results/${study_name}_04.tgz" ]; then
        echo "hase tar results present"
    else
        echo "hase tar results absent. Please re-run"
    fi

}


check_results_07 () {

    #chunk_count=`cat ${section_07_dir}/tabfile.info1 | wc -l`
    #BFfile=`ls ${section_07_dir}/vQTL_BF_*besd | wc -l`
    #svlmfile=`ls ${section_07_dir}/vQTL_BF_*besd | wc -l`
    #drmfile=`ls ${section_07_dir}/vQTL_BF_*besd | wc -l`

    #if [ $chunk_count = $BFfile ]; then
    #    echo "vQTL detection results with BF method present"
    #fi

    #if [ $chunk_count = $svlmfile ]; then
    #    echo "vQTL detection results with svlm method present"
    #fi 

    #if [ $chunk_count = $drmfile ]; then
    #    echo "vQTL detection results with drm method present"
    #fi
    
    tarfile=`ls ${home_directory}/results/${study_name}_07_chr*.tgz | wc -l`
    if [ $tarfile = 22 ]; then
        echo "vQTL tar results present"
    else
        echo "vQTL tar results are missing. Please re-run"
    fi
}


check_results_08 () {

	if [ -f "${section_08_dir}/badInversions.txt" ]; then
		echo "Bad inversions file present"
	else
		echo "Problem: Bad inversions file is absent"
		exit 1
	fi

	if [ -f "${section_08_dir}/inversionsSummary.txt" ]; then
		echo "Inversion frequency file present"
	else
		echo "Problem: Inversion frequency file is absent"
		exit 1
	fi

	if [ -f "${section_08_dir}/invmeqtl.Rdata" ]; then
		echo "inversionmeQTL statistics file present"
	else
		echo "Problem: inversionmeQTL statistics file is absent"
		exit 1
	fi

}



check_results_09 () {


  vect_PRS=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"_" '{print $2}' |tr "\n" " ")
  vect_PRS_array=($vect_PRS)

  n=$((${#vect_PRS_array[*]}-1))

  for ((k=0;k<=$n;k++))
  do
    PRS=${vect_PRS_array[$k]}
    PRS_dir=${section_09_dir}/${PRS}
   
	  if [ -f "${section_09_dir}/${PRS}/${study_name}_PRS_${PRS}_hist.pdf" ]; then
		  echo "Histogram for $PRS PRS present"
	  else
		  echo "Problem: histogram for $PRS PRS is absent"
		  exit 1
	  fi

	  if [ -f "${section_09_dir}/${PRS}/${study_name}_PRS_${PRS}_cell_counts_plots.pdf" ]; then
		  echo "Cell count PRS correlation plot for $PRS PRS present"
	  else
		  echo "Problem: cell count correlation plot for $PRS PRS is absent"
		  exit 1
	  fi

	  if [ -f "${section_09_dir}/${PRS}/${study_name}_PRS_${PRS}_EWAS_results.RData" ]; then
		  echo "EWAS results file for $PRS PRS present"
	  else
		  echo "Problem: EWAS results file for $PRS PRS is absent"
		  exit 1
	  fi

	  if [ -f "${section_09_dir}/${PRS}/${study_name}_PRS_${PRS}_EWAS.ewas.report.html" ]; then
		  echo "EWAS report file for $PRS PRS present"
	  else
		  echo "Problem: EWAS report file for $PRS PRS is absent"
		  exit 1
	  fi

		 if [ -f "${section_09_dir}/${PRS}/${study_name}_PRS_${PRS}_EWAS_qqplot.nocovs.pdf" ]; then
		
      echo "QQplots for $PRS PRS EWAS present"
	  else
		  echo "Problem: not all QQplots for $PRS PRS EWAS present"
		  exit 1
	  fi

  done

}

check_results_14 () {

	if [ -f "${section_14_dir}/nc886_scatter.jpeg" ]; then
		echo "scatterplot present"
	else
		echo "Problem: Scatterplot is absent"
		exit 1
	fi

	if [ -f "${section_14_dir}/nc886_frequency.txt" ]; then
		echo "nc866 frequency file present"
	else
		echo "Problem: nc866 frequency file is absent"
		exit 1
	fi
 
	if [ -f "${section_14_dir}/nc886_groups.txt" ]; then
		rm ${section_14_dir}/nc886_groups.txt
  
	fi
	}

ld_manifest_chunks_15 () {
	local outdir="$1" chr="$2"
	python - "${outdir}" "${chr}" <<'PY'
import json
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
chrom = sys.argv[2]
manifest = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
actual_filter = manifest.get("variant_index", {}).get("chromosome_filter")
if str(actual_filter) != str(chrom):
    raise SystemExit(
        f"manifest chromosome_filter is {actual_filter!r}, expected {chrom!r}"
    )
chrom_meta = manifest.get("A_blocks", {}).get("chromosomes", {}).get(str(chrom))
if not chrom_meta:
    raise SystemExit(f"manifest lacks A_blocks metadata for chr{chrom}")
chunks = [chunk.get("name") for chunk in chrom_meta.get("chunks", [])]
if any(not name for name in chunks):
    raise SystemExit(f"manifest has an unnamed A-block chunk for chr{chrom}")
if int(chrom_meta.get("n_chunks", -1)) != len(chunks):
    raise SystemExit(f"manifest n_chunks does not match chunk list for chr{chrom}")
if not chunks:
    raise SystemExit(f"manifest lists no A-block chunks for chr{chrom}")
for name in chunks:
    print(name)
PY
}

ld_verify_manifest_checksum_15 () {
	local outdir="$1"
	PYTHONPATH="${scripts_directory}/resources/genetics:${PYTHONPATH:-}" python - "${outdir}" <<'PY'
import sys
from pathlib import Path
import ld_checksums

outdir = Path(sys.argv[1])
document = ld_checksums.read_cohort_checksums(outdir)
expected = (document.get("files") or {}).get("manifest.json")
if not expected:
    raise SystemExit("checksums.json does not contain manifest.json")
actual = ld_checksums.hash_file(outdir / "manifest.json")
if actual != expected:
    raise SystemExit(
        f"manifest.json checksum mismatch: expected {expected}, got {actual}"
    )
PY
}

check_section_15_central_results () {

	if [ -d "${ld_precursor_dir}" ]; then
		if [ -f "${ld_precursor_dir}/precursor_manifest.json" ]; then
			echo "LD precursor manifest present"
		else
			echo "Problem: LD precursor manifest is absent"
			exit 1
		fi
	fi

	if [ -d "${ld_panel_dir}" ]; then
		latest_panel=""
		latest_panel_n=0
		panel_has_entries=0
		shopt -s nullglob
		for panel_entry in "${ld_panel_dir}"/*; do
			panel_has_entries=1
		done
		for panel in "${ld_panel_dir}"/panel_v*; do
			panel_n="${panel##*panel_v}"
			if [[ "${panel_n}" =~ ^[0-9]+$ ]] && [ "${panel_n}" -gt "${latest_panel_n}" ]; then
				latest_panel="${panel}"
				latest_panel_n="${panel_n}"
			fi
		done
		shopt -u nullglob
		if [ -n "${latest_panel}" ] && [ -f "${latest_panel}/pooled_manifest.json" ]; then
				echo "LD pooled panel manifest present in ${latest_panel}"
		elif [ -f "${ld_panel_dir}/pooled_manifest.json" ]; then
			echo "LD pooled panel manifest present"
		elif [ "${panel_has_entries}" -eq 1 ]; then
			echo "Problem: LD pooled panel manifest is absent"
			exit 1
		fi
	fi

}

check_section_15_upload_dir () {
	local path base ok=1 section_15_upload_dir="${section_15_dir}/upload"
	if [ ! -d "${section_15_upload_dir}" ]; then
		echo "Problem: LD upload staging directory is absent: ${section_15_upload_dir}"
		exit 1
	fi
	shopt -s nullglob dotglob
	for path in "${section_15_upload_dir}"/*; do
		base="$(basename "${path}")"
		case "${base}" in
			*.tgz.aes|*.md5sum|.uploaded_*) ;;
			*)
				echo "Problem: unexpected raw or unsupported file in LD upload directory: ${path}"
				exit 1
				;;
		esac
	done
	shopt -u dotglob
	for path in "${section_15_upload_dir}"/*.tgz.aes; do
		ok=0
		if [ ! -f "${path%.tgz.aes}.md5sum" ]; then
			echo "Problem: missing md5sum for ${path}"
			exit 1
		fi
	done
	shopt -u nullglob
	if [ "${ok}" -ne 0 ]; then
		echo "Problem: no LD encrypted archives found in ${section_15_upload_dir}"
		exit 1
	fi
}

check_results_15 () {
	local chunks targets
	if ! targets="$(ld_target_chromosomes_15)"; then
		exit 1
	fi

	check_section_15_upload_dir
	section_15_upload_dir="${section_15_dir}/upload"

	while IFS= read -r chr; do
		outdir="$(ld_resolve_chromosome_dir_15 "${ld_prepare_dir}" "${chr}")"
		if ! ld_verify_manifest_checksum_15 "${outdir}"; then
			echo "Problem: LD chr${chr} manifest does not match checksums.json"
			exit 1
		fi
		if [ ! -f "${outdir}/.packaged" ]; then
			echo "Problem: LD chr${chr} packaged sentinel is absent"
			exit 1
		fi
		for f in manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt; do
			if [ ! -f "${outdir}/${f}" ]; then
				echo "Problem: LD chr${chr} ${f} is absent"
				exit 1
			fi
		done
	for path in "${section_15_upload_dir}"/*.md5sum; do
		if [ ! -f "${path%.md5sum}.tgz.aes" ]; then
			echo "Problem: missing encrypted archive for ${path}"
			exit 1
		fi
	done
		if [ -d "${outdir}/A_blocks" ]; then
			echo "Problem: LD chr${chr} A_blocks still exists after packaging"
			exit 1
		fi
		scaffold="${study_name}_chr${chr}_15_scaffold"
		if ! ld_verify_archive "${section_15_upload_dir}" "${scaffold}"; then
			echo "Problem: LD chr${chr} scaffold upload artefacts are absent"
			exit 1
		fi
		chunk_count=0
		if ! chunks="$(ld_manifest_chunks_15 "${outdir}" "${chr}")"; then
			echo "Problem: LD chr${chr} chunk manifest is invalid"
			exit 1
		fi
		while IFS= read -r chunk; do
			chunk_base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
			if ! ld_verify_archive "${section_15_upload_dir}" "${chunk_base}"; then
				echo "Problem: LD chr${chr} chunk upload artefacts are absent for ${chunk}"
				exit 1
			fi
			chunk_count=$((chunk_count + 1))
		done <<< "${chunks}"
		if [ "${chunk_count}" -eq 0 ]; then
			echo "Problem: LD chr${chr} has no chunk upload artefacts"
			exit 1
		fi
		echo "LD chr${chr} packaged artefacts present"
	done <<< "${targets}"

}
