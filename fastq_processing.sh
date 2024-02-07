#!/bin/bash

# Input użytkownika, którym jest nunmer projektu
read -p 'Enter the number of BioProject: ' no_project

if [[ -n "$no_project" ]]; then

        # Krok 1 - pobieranie informacji o projeckie podanym przez użytkownika
        esearch -db sra -query "$no_project" | efetch -format runinfo -mode xml | xtract -pattern SraRunInfo -element R>

        ## obróbka

        tr "\t" "\n" < runinfo.txt > SRR_list.txt

        ## Przygotowanie folderów do podziału na pair_end and single_end

        mkdir pair_end single_end

        # Krok 2
        read -p 'Enter the number of readings for the fastq-dump programme: ' no_reads

        while IFS= read -r AR; do
                echo "*** Pobieranie: ${AR}";
                fastq-dump -X $no_reads ${AR} --split-files

                if [ -f ${AR}_1.fastq ] && [ -f ${AR}_2.fastq ]; then
                    mv ${AR}_1.fastq pair_end
                    mv ${AR}_2.fastq pair_end
                else
                    mv ${AR}_1.fastq single_end
                fi

        done < SRR_list.txt

        # Krok 3 kontrola jakości

        mkdir post_trimm pre_trimm

        fastqc pair_end/* -o pre_trimm

        fastqc single_end/* -o pre_trimm

        # Krok 4 trimming

        mkdir trimmed_pair_end trimmed_single_end

        # Trimming of pair end records

        seq_1_list=(./pair_end/*_1*)
        seq_2_list=(./pair_end/*_2*)

        for i in "${!seq_1_list[@]}"; do
                filename1=$(basename -- "${seq_1_list[i]}");
                filename1="${filename1%.*}";
                filename2=$(basename -- "${seq_2_list[i]}");
                filename2="${filename2%.*}";

                echo Trimming of $filename1 and $filename2 ...;
               java -jar /bioapp/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${seq_1_list[i]} ${seq_2_list[i]} trimmed_$filename1.fq unpaired_$filename1.fq trimmed_$filename2.fq unpaired_$filena>

                mv trimmed_*.fq ./trimmed_pair_end;
                rm ./unpaired*;
        done

        # Trimming of single end records

        single_list=./single_end/*.fastq

        for file in $single_list; do
                filename=$(basename -- "${file}");
                filename="${filename%.*}";

                echo Trimming of $filename...;
                java -jar /bioapp/Trimmomatic-0.39/trimmomatic-0.39.jar SE $file trimmed_$filename.fq SLIDINGWINDOW:4:30

                mv trimmed_*.fq ./trimmed_single_end;
        done

        # Krok 5 - quality control of trimmed files

        fastqc ./trimmed_pair_end/* -o post_trimm

        fastqc ./trimmed_single_end/* -o post_trimm

        # Krok 6 - mapowanie z genomem hg19

        # pobieranie pliku z genomem referencyjnym hg19
        #wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz

        # rozpakowywanie pliku hg19
        #gzip -dk hg19.fa.gz

        # Tworzenie indeksu dla programu hisat2 na podsatwie genomu referencyjnego
        #hisat2-build -p 4 -f hg19.fa hg19

        #mkdir ./refs

        #mv hg19.*.ht2 hg19.fa.gz hg.19fa ./refs

        # mapowanie sekwencji pair end do genomu referencyjnego

        mkdir sam_files

        seq_1_list=($PWD/trimmed_pair_end/*_1.fq)
        seq_2_list=($PWD/trimmed_pair_end/*_2.fq)

        for i in "${!seq_1_list[@]}"; do
                filename=$(basename -- "${seq_1_list[i]}");
                filename="${filename%_*}";
                echo $filename;

                hisat2 -x ../refs/hg19 -1 ${seq_1_list[i]} -2 ${seq_2_list[i]} -S $filename.sam;
        done

        mv *.sam ./sam_files

        # mapowanie sekwencji single end do genomu referencyjnego
        single_seq=$PWD/trimmed_single_end/*.fq

        for file in $single_seq; do
                filename=$(basename -- "${file}");
                filename="${filename%_*}";

                hisat2 -q -x ../refs/hg19 -U $file -S $filename.sam;
        done

        mv *.sam ./sam_files

        # Krok 7 - konwersja plików sam do bam

        mkdir bam_files

        sams=$PWD/sam_files/*.sam

        for fsam in $sams; do
                filename=$(basename -- "${fsam}");
                filename="${filename%.*}";

                samtools view -Sb -@ 4 $fsam > $filename.bam | mv $filename.bam ./bam_files;
        done

        # Krok 8 - zliczanie odczytów

        # Pobieranie pliku gtf dla odpowiedniego genomu
        # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz

        # Rozpakowanie i zmiana nazwy
        # gzip -dk hg19.refGene.gtf.gz | mv hg19.reFGene.gtf hg19.gtf | mv ./hg19 ../refs

        featureCounts -a ../refs/hg19.gtf -g gene_name -o counts.txt ./bam_files/*.bam

else
        echo "Nie podano numeru projektu";
fi;
