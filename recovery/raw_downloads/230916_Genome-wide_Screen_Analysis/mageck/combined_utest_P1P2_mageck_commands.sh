## 09/20/23
## run in local computer
cd /Users/dulguun/oak/Users/dulguun/230916_Genome-wide_Screen_Analysis
conda activate mageckenv
mageck test -k from_chad/combined_utest_P1P2_mageck_input.txt -t Pos -c Neg --control-sgrna from_chad/mageck/non_targeting_sgRNAs.txt -n Pos_vs_Neg
mageck test -k from_chad/combined_utest_P1P2_mageck_input.txt -t EC_total -c D0 --control-sgrna from_chad/mageck/non_targeting_sgRNAs.txt -n EC_total_vs_D0


## copied from Chad's previous run
mageck test -k combined_utest_P1P2_mageck_input.txt -t D0 -c Pos --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/D0_vs_Pos/D0_vs_Pos
mageck test -k combined_utest_P1P2_mageck_input.txt -t D0 -c Neg --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/D0_vs_Neg/D0_vs_Neg
mageck test -k combined_utest_P1P2_mageck_input.txt -t D0 -c EC_total --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/D0_vs_EC_total/D0_vs_EC_total
mageck test -k combined_utest_P1P2_mageck_input.txt -t D0 -c Skew --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/D0_vs_Skew/D0_vs_Skew
mageck test -k combined_utest_P1P2_mageck_input.txt -t Pos -c Neg --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/Pos_vs_Neg/Pos_vs_Neg
mageck test -k combined_utest_P1P2_mageck_input.txt -t Pos -c EC_total --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/Pos_vs_EC_total/Pos_vs_EC_total
mageck test -k combined_utest_P1P2_mageck_input.txt -t Pos -c Skew --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/Pos_vs_Skew/Pos_vs_Skew
mageck test -k combined_utest_P1P2_mageck_input.txt -t Neg -c EC_total --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/Neg_vs_EC_total/Neg_vs_EC_total
mageck test -k combined_utest_P1P2_mageck_input.txt -t Neg -c Skew --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/Neg_vs_Skew/Neg_vs_Skew
mageck test -k combined_utest_P1P2_mageck_input.txt -t EC_total -c Skew --control-sgrna Data/guideLibs/non_targeting_sgRNAs.txt -n combined_utest_P1P2/mageck/EC_total_vs_Skew/EC_total_vs_Skew
