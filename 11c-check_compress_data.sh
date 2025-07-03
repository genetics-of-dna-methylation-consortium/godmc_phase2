#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_11_dir}/logs_c
touch ${section_11c_logfile}
exec &> >(tee ${section_11c_logfile})
print_version

#### fucntions for version comparsion ############## 
source resources/logs/check_logs.sh


#### fucntions for checking if the files has been updated ############## 
check_file_updated () {
	file=$1
	script_name=$2
	if find "$file" -newermt "2025-07-01" | grep -q .; then
    	echo "File $file was modified after 01 July."
	else
		echo "Problem: File $file was NOT modified after 01 July."
		echo "Please rerun the corresponding scripts $script_name after git pull."
		exit 1
	fi
}


#### fucntions for checking log files ############## 

check_logs () {
	file=$1
	script_name=$2
	# compare_version $3
	if grep -i -q "success" "$file"; then
		echo "$script_name completed successfully."
	else
		echo "Problem: $script_name did not complete successfully"
		exit 1
	fi
}

#### fucntions for checking output files ############## 

check_results_10a () {

	if [ -f "${section_03_dir}/age_prediction.pdf" ]; then
		echo "Plot of age predictions is present"
	else
		echo "Problem: the plot of age predictions is absent"
		exit 1
	fi

	check_file_updated "${section_03_dir}/age_prediction.pdf" "10a"
	
	if [ -f "${section_03_dir}/age_prediction_correlation.png" ]; then
		echo "Plot of correleation on age predictions is present"
	else
		echo "Problem: the plot of age predictions is absent"
		exit 1
	fi

	check_file_updated "${section_03_dir}/age_prediction_correlation.png" "03a"

	if [ -f "${section_03_dir}/age_prediction_stats.csv" ]; then
		echo "Statistic table of age predictions is present"
	else
		echo "Problem: Statistic table of of age predictions is absent"
		exit 1
	fi

	check_file_updated "${section_03_dir}/age_prediction_stats.csv" "10a"

	if [ -f "${section_03_dir}/age_prediction_stats_models.csv" ]; then
		echo "Statistic table of age predictions of models is present"
	else
		echo "Problem: Statistic table of of age predictions of models is absent"
		exit 1
	fi

	check_file_updated "${section_03_dir}/age_prediction_stats_models.csv" "10a"

	if [ -f "${section_03_dir}/age_prediction_stats_corrsd.csv" ]; then
		echo "Correleation table of age predictions is present"
	else
		echo "Problem: Correleation table of of age predictions is absent"
		exit 1
	fi

	check_file_updated "${section_03_dir}/age_prediction_stats_corrsd.csv" "10a"

	GWASresults=$(find ${section_10_dir} -type f -name "*.fastGWA" | wc -l)
	if [ $GWASresults -gt 10 ]; then
		echo "The GWAS results are present"
	else
		echo "Problem: GWAS has not been performed successfully"
	    exit 1
	fi

	for file in $(find ${section_10_dir} -type f -name "*.fastGWA");
	do
		check_file_updated "$file" "10a"
	done

	MANplot=$(find ${section_10_dir} -type f -name "*_manhattan.pdf" | wc -l)
	if [ $MANplot -gt 4 ]; then
		echo "Manhattan plots of age accelerations are present"
	else
		echo "Problem: Manhattan plots of age accelerations are absent"
		exit 1
	fi

	for file in $(find ${section_10_dir} -type f -name "*_manhattan.pdf");
	do
		check_file_updated "$file" "10a"
	done

	QQplot=$(find ${section_10_dir} -type f -name "*_qqplot.jpeg" | wc -l)
	if [ $QQplot -gt 4 ]; then
		echo "QQ plots of age accelerations are present"
	else
		echo "Problem: QQ plots of age accelerations are absent"
		exit 1
	fi

	for file in $(find ${section_10_dir} -type f -name "*_qqplot.jpeg");
	do
		check_file_updated "$file" "10a"
	done

	
}

check_results_10b () {

	HERresults=$(find ${section_10_dir} -type f -name "heritability_*.hsq" | wc -l)
	if [ $HERresults -ge 2 ]; then
		echo "Hertabililty results of age predictions are present"
	else
		echo "Problem: Heritability calculation for age acceleration hasn't been performed successfully"
		exit 1
	fi

	for file in $(find ${section_10_dir} -type f -name "heritability_*.hsq");
	do
		check_file_updated "$file" "10b"
	done

}

check_results_11a () {

	if [ -f "${smoking_pred_plot}" ]; then
		echo "Smoking prediction plot is present"
	else
		echo "Problem: Smoking prediction plot is not present"
		exit 1
	fi

	check_file_updated "${smoking_pred_plot}" "11a"

	if [ -f "${home_directory}/results/11/smoking_stats.csv" ]; then
		echo "Statistic table of preicted smoking is present"
	else
		echo "Problem: Statistic table of preicted smoking is not present"
		exit 1
	fi

	check_file_updated "${home_directory}/results/11/smoking_stats.csv" "11a"

	if [ -f "${section_11_dir}/gwas_smoking.fastGWA" ] && [ -f "${section_11_dir}/gwas_smoking_PCA.fastGWA"  ]; then
		echo "The GWAS result of smoking is present"
	else
		echo "Problem: GWAS has not been performed successfully"
		exit 1
	fi
	
	check_file_updated "${section_11_dir}/gwas_smoking.fastGWA" "11a"
	check_file_updated "${section_11_dir}/gwas_smoking_PCA.fastGWA" "11a"

	if [ -f "${section_11_dir}/gwas_smoking_manhattan.pdf" ] && [ -f "${section_11_dir}/gwas_smoking_PCA_manhattan.pdf" ]; then
		echo "Manhattan plot of smoking is present"
	else
		echo "Problem: Manhattan plot of smoking is absent"
		exit 1
	fi

	check_file_updated "${section_11_dir}/gwas_smoking_manhattan.pdf" "11a"
	check_file_updated "${section_11_dir}/gwas_smoking_PCA_manhattan.pdf" "11a"
	
	
	if [ -f "${section_11_dir}/gwas_smoking_qqplot.jpeg" ] && [ -f "${section_11_dir}/gwas_smoking_PCA_qqplot.jpeg" ]; then
		echo "QQ plot of smoking is present"
	else
		echo "Problem: QQ plot of smoking is absent"
		exit 1
	fi

	check_file_updated "${section_11_dir}/gwas_smoking_qqplot.jpeg" "11a"
	check_file_updated "${section_11_dir}/gwas_smoking_PCA_qqplot.jpeg" "11a"

	
}

check_results_11b () {	
	if [ -f "${section_11_dir}/heritability_smoking.hsq" ]; then
		echo "Hertabililty result of smoking is present"
	else
		echo "Problem: Heritability calculation for smoking hasn't been performed successfully"
		exit 1
	fi

	check_file_updated "${section_11_dir}/heritability_smoking.hsq" "11b"
	
}


### check the log files and the output files for 03a, 10 and 11 ############## 
script_names=( "10a-gwas_aar.sh" "10b-gwas_aar.sh" "11a-gwas_smoking.sh" "11b-gwas_smoking.sh" )
logfile_names=( "${section_10a_logfile}" "${section_10b_logfile}" "${section_11a_logfile}" "${section_11b_logfile}" )
script_number=( "10a" "10b" "11a" "11b" )


for index in "${!script_names[@]}"
do
	script=${script_names[$index]}
	log=${logfile_names[$index]}
	num_script=${script_number[$index]}

	echo ""
  	echo "Checking log files for $script"

	check_logs ${log} ${script} ${num_script}

  	echo "Checking results files for $script"
	
	check_results_${num_script}

	echo "Section $num_script has been successfully completed!"
done


#### Compress the output files from 03a, 10a, 10b, 11a and 11b ############## 
echo ""
echo "Compressing the outputs from 03a, 10 and 11"
cd ${home_directory}

tar -zcf results/AgeSmokGWAS2025_${study_name}.tgz results/03/logs_a/log.txt results/03/age_prediction.pdf results/03/age_prediction_correlation.png results/03/age_prediction_stats.csv results/03/age_prediction_stats_corrsd.csv results/03/smoking_prediction.pdf results/03/cohort_descriptives_commonids.RData results/10 results/11 
tar -zcf results/AgeSmokGWAStest_${study_name}.tgz results/03/logs_a/log.txt results/03/age_prediction* results/03/smoking_prediction.pdf results/03/cohort_descriptives_commonids.RData results/10 results/11 

echo "Successfully created results archives ${home_directory}/results/11/AgeSmokGWAS2025_${study_name}.tgz"


# Generating md5 checksum for verify the data intensity
cd ${home_directory}/results || exit 1
md5sum AgeSmokGWAS2025_${study_name}.tgz > AgeSmokGWAS2025_${study_name}.tgz.md5sum
md5sum -c AgeSmokGWAS2025_${study_name}.tgz.md5sum
# encryption 
gpg --output AgeSmokGWAS2025_${study_name}.tgz.gpg --symmetric --cipher-algo AES256 AgeSmokGWAS2025_${study_name}.tgz
echo ""
#echo "Please download the follwing files to your own local machine and upload via the link of https://drive.google.com/drive/folders/19T0aSzh7xX6rX17pe6HBNImMYjoZEGSW?usp=sharing."
echo "!!!!!!Please noted: the link of GoogleDrive has been updated and the old one will not work any more!!!!!!"
echo "!!!!!!Please noted: Please upload the file named with 2025 which is the updated folder!!!!!!"
echo "Please download the following files to your own local machine and upload via the link of https://drive.google.com/drive/folders/1CvrU4qDNJSS2J8a_MtlGm33UCDzVby5Y?usp=sharing."
echo "1. " ${home_directory}/results/AgeSmokGWAS2025_${study_name}.tgz.md5sum
echo "2. " ${home_directory}/results/AgeSmokGWAS2025_${study_name}.tgz.gpg
echo "Please share encryption passphrase to the developers by emailing the developers s.w.wang@exeter.ac.uk mentioned in the wiki."
echo "Thank you very much. Have a nice day!"
# decompressing archives
# gpg -d -o AgeSmokGWAS.tgz AgeSmokGWAS.tgz.gpg
# tar -xzvf AgeSmokGWAS.tgz


