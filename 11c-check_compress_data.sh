#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_11_dir}/logs_c
touch ${section_11c_logfile}
exec &> >(tee ${section_11c_logfile})
print_version

#### fucntions for version comparsion ############## 
source resources/logs/check_logs.sh

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
check_results_03a () {
	if [ -f "${section_03_dir}/age_prediction.pdf" ]; then
		echo "Plot of age predictions is present"
	else
		echo "Problem: the plot of age predictions is absent"
		exit 1
	fi

	
	if [ -f "${section_03_dir}/age_prediction_correlation.png" ]; then
		echo "Plot of correleation on age predictions is present"
	else
		echo "Problem: the plot of age predictions is absent"
		exit 1
	fi

	if [ -f "${section_03_dir}/age_prediction_stats.csv" ]; then
		echo "Statistic table of age predictions is present"
	else
		echo "Problem: Statistic table of of age predictions is absent"
		exit 1
	fi

	if [ -f "${section_03_dir}/age_prediction_stats_corrsd.csv" ]; then
		echo "Correleation table of age predictions is present"
	else
		echo "Problem: Correleation table of of age predictions is absent"
		exit 1
	fi

}

check_results_10a () {

	GWASresults=$(find ${section_10_dir} -type f -name "*.fastGWA" | wc -l)
	if [ $GWASresults -gt 2 ]; then
		echo "The GWAS results are present"
	else
		echo "Problem: GWAS has not been performed successfully"
	    exit 1
	fi

	MANplot=$(find ${section_10_dir} -type f -name "*_manhattan.pdf" | wc -l)
	if [ $MANplot -gt 2 ]; then
		echo "Manhattan plots of age accelerations are present"
	else
		echo "Problem: Manhattan plots of age accelerations are absent"
		exit 1
	fi

	QQplot=$(find ${section_10_dir} -type f -name "*_qqplot.png" | wc -l)
	if [ $QQplot -gt 2 ]; then
		echo "QQ plots of age accelerations are present"
	else
		echo "Problem: QQ plots of age accelerations are absent"
		exit 1
	fi

}

check_results_10b () {

	HERresults=$(find ${section_10_dir} -type f -name "heritability_*.hsq" | wc -l)
	if [ $HERresults -ge 2 ]; then
		echo "Hertabililty results of age predictions are present"
	else
		echo "Problem: Heritability calculation for age acceleration hasn't been performed successfully"
		exit 1
	fi

}

check_results_11a () {
	if [ -f "${section_11_dir}/gwas_smoking.fastGWA" ]; then
		echo "The GWAS result of smoking is present"
	else
		echo "Problem: GWAS has not been performed successfully"
		exit 1
	fi
	
	if [ -f "${section_11_dir}/gwas_smoking_manhattan.pdf" ]; then
		echo "Manhattan plot of smoking is present"
	else
		echo "Problem: Manhattan plot of smoking is absent"
		exit 1
	fi
	
	
	if [ -f "${section_11_dir}/gwas_smoking_qqplot.png" ]; then
		echo "QQ plot of smoking is present"
	else
		echo "Problem: QQ plot of smoking is absent"
		exit 1
	fi
	
}

check_results_11b () {	
	if [ -f "${section_11_dir}/heritability_smoking.hsq" ]; then
		echo "Hertabililty result of smoking is present"
	else
		echo "Problem: Heritability calculation for smoking hasn't been performed successfully"
		exit 1
	fi

	
}


### check the log files and the output files for 03a, 10 and 11 ############## 
script_names=( "03a-methylation_variables.sh" "10a-gwas_aar.sh" "10b-gwas_aar.sh" "11a-gwas_smoking.sh" "11b-gwas_smoking.sh" )
logfile_names=( "${section_03a_logfile}"  "${section_10a_logfile}" "${section_10b_logfile}" "${section_11a_logfile}" "${section_11b_logfile}" )
script_number=( "03a" "10a" "10b" "11a" "11b" )


for index in "${!script_names[@]}"
do
	script=${script_names[$index]}
	log=${logfile_names[$index]}
	num_script=${script_number[index]}

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
tar -zcf results/AgeSmokGWAS_${study_name}.tgz results/03/age_prediction.pdf results/03/age_prediction_correlation.png results/03/age_prediction_stats.csv results/03/age_prediction_stats_corrsd.csv results/10 results/11 
					
echo "Successfully created results archives ${home_directory}/results/11/AgeSmokGWAS_${study_name}.tgz"


# Generating md5 checksum for verify the data intensity
cd ./results
md5sum AgeSmokGWAS_${study_name}.tgz > AgeSmokGWAS_${study_name}.tgz.md5sum
md5sum -c AgeSmokGWAS_${study_name}.tgz.md5sum
# encryption 
gpg --output AgeSmokGWAS_${study_name}.tgz.gpg --symmetric --cipher-algo AES256 AgeSmokGWAS_${study_name}.tgz
echo ""
echo "Please download the follwing files to your own local machine and upload via the link of https://drive.google.com/drive/folders/19T0aSzh7xX6rX17pe6HBNImMYjoZEGSW?usp=sharing."
echo "1. " ${home_directory}/results/AgeSmokGWAS_${study_name}.tgz.md5sum
echo "2. " ${home_directory}/results/AgeSmokGWAS_${study_name}.tgz.gpg
echo "Please share encryption passphrase to the developers by emailing the developers mentioned in the wiki."
echo "Thank you very much. Have a nice day!"
# decompressing archives
# gpg -d -o AgeSmokGWAS.tgz AgeSmokGWAS.tgz.gpg
# tar -xzvf AgeSmokGWAS.tgz


