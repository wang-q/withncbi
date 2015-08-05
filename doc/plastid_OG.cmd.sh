mkdir -p ~/data/organelle/plastid_OG
cd ~/data/organelle/plastid_OG

#----------------------------------------------------------#
# Green alga TODO
#----------------------------------------------------------#
# Chlorella
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Chlorella \
    -t Chlore_sorokiniana \
    -q Chlore_vulgaris \
    -q Chlore_mirabilis \
    -q Chlore_variabilis

# Koliella
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Koliella \
    -t Ko_longiseta \
    -q Ko_corcontica

#----------------------------------------------------------#
# Angiosperm
#----------------------------------------------------------#
#----------------------------#
# Araliaceae 五加科
#----------------------------#
# Dendropanax
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Dendropanax \
    -t Dendrop_morbifer \
    -q Dendrop_dentiger \
    -q Aral_undulata \
    -o Aral_undulata

# Panax
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Panax \
    -t Panax_ginseng \
    -q Panax_notoginseng \
    -q Panax_quinquefolius \
    -q Aral_undulata \
    -o Aral_undulata

#----------------------------#
# Orchidaceae 兰科
#----------------------------#
# Corallorhiza
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Corallorhiza \
    -t Cora_trifida \
    -q Cora_bulbosa \
    -q Cora_mertensiana \
    -q Cora_odontorhiza \
    -q Cora_wisteriana \
    -q Cora_macrantha \
    -q Catt_crispata  \
    -o Catt_crispata

# Cymbidium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cymbidium \
    -t Cym_aloifolium \
    -q Cym_sinense \
    -q Cym_tracyanum \
    -q Cym_tortisepalum \
    -q Cym_mannii \
    -q Catt_crispata \
    -o Catt_crispata

# Cypripedium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cypripedium \
    -t Cyp_formosanum \
    -q Cyp_macranthos \
    -q Cyp_japonicum \
    -q Catt_crispata \
    -o Catt_crispata

# Epipogium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Epipogium \
    -t Epip_aphyllum \
    -q Epip_roseum \
    -q Catt_crispata \
    -o Catt_crispata

# Masdevallia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Masdevallia \
    -t Mas_coccinea \
    -q Mas_picturata \
    -q Catt_crispata \
    -o Catt_crispata

# Paphiopedilum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Paphiopedilum \
    -t Pap_armeniacum \
    -q Pap_niveum \
    -q Catt_crispata \
    -o Catt_crispata

# Phalaenopsis
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Phalaenopsis \
    -t Phalae_equestris \
    -q Phalae_hybrid_cultivar \
    -q Phalae_aphrodite_subsp_formosana \
    -q Catt_crispata \
    -o Catt_crispata

#----------------------------#
# Asteraceae 菊科
#----------------------------#
# Artemisia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Artemisia \
    -t Art_frigida \
    -q Art_montana \
    -q Gui_abyssinica \
    -o Gui_abyssinica

# Chrysanthemum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Chrysanthemum \
    -t Chrysa_x_morifolium \
    -q Chrysa_indicum \
    -q Gui_abyssinica \
    -o Gui_abyssinica

# Helianthus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Helianthus \
    -t Helia_annuus \
    -q Helia_tuberosus \
    -q Helia_decapetalus \
    -q Helia_divaricatus \
    -q Helia_giganteus \
    -q Helia_grosseserratus \
    -q Helia_hirsutus \
    -q Helia_maximiliani \
    -q Helia_strumosus \
    -q Gui_abyssinica \
    -o Gui_abyssinica

#----------------------------#
# Brassicaceae 十字花科
#----------------------------#
# Aethionema
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Aethionema \
    -t Aet_grandiflorum \
    -q Aet_cordifolium \
    -q Lep_virginicum \
    -o Lep_virginicum

# Arabis
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Arabis \
    -t Arabis_alpina \
    -q Arabis_hirsuta \
    -q Lep_virginicum \
    -o Lep_virginicum

# Brassica
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Brassica \
    -t Brassi_napus \
    -q Brassi_rapa_subsp_pekinensis \
    -q Lep_virginicum \
    -o Lep_virginicum

# Capsella
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Capsella \
    -t Capse_bursa_pastoris \
    -q Capse_rubella \
    -q Lep_virginicum \
    -o Lep_virginicum

# Cardamine
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cardamine \
    -t Car_impatiens \
    -q Car_resedifolia \
    -q Lep_virginicum \
    -o Lep_virginicum

# Pachycladon
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pachycladon \
    -t Pachyc_cheesemanii \
    -q Pachyc_enysii \
    -q Lep_virginicum \
    -o Lep_virginicum

#----------------------------#
# Amaranthaceae 苋科
#----------------------------#
# Salicornia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Salicornia \
    -t Salic_bigelovii \
    -q Salic_brachiata \
    -q Salic_europaea \
    -q Spi_oleracea \
    -o Spi_oleracea

#----------------------------#
# Caryophyllaceae 石竹科
#----------------------------#
# Silene
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Silene \
    -t Si_latifolia \
    -q Si_chalcedonica \
    -q Si_conica \
    -q Si_noctiflora \
    -q Si_vulgaris \
    -q Si_paradoxa \
    -q Si_conoidea \
    -q Agroste_githago \
    -o Agroste_githago

#----------------------------#
# Fabaceae 豆科
#----------------------------#
# Glycine
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Glycine \
    -t Glyci_max \
    -q Glyci_soja \
    -q Glyci_tomentella \
    -q Glyci_cyrtoloba \
    -q Glyci_falcata \
    -q Glyci_canescens \
    -q Glyci_dolichocarpa \
    -q Glyci_stenophita \
    -q Glyci_syndetika \
    -q Phas_vulgaris \
    -o Phas_vulgaris

# Lathyrus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Lathyrus \
    -t Lat_clymenum \
    -q Lat_odoratus \
    -q Lat_sativus \
    -q Lat_tingitanus \
    -q Lat_japonicus \
    -q Lat_davidii \
    -q Lat_littoralis \
    -q Lat_palustris \
    -q Lat_pubescens \
    -q Lat_graminifolius \
    -q Lat_ochroleucus \
    -q Lat_inconspicuus \
    -q Lat_venosus \
    -q Cic_arietinum \
    -o Cic_arietinum

# Lupinus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Lupinus \
    -t Lu_albus \
    -q Lu_luteus \
    -q Glycy_glabra \
    -o Glycy_glabra

# Medicago
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Medicago \
    -t Med_truncatula \
    -q Med_hybrida \
    -q Med_papillosa \
    -q Cic_arietinum \
    -o Cic_arietinum

# Trifolium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Trifolium \
    -t Trif_repens \
    -q Trif_subterraneum \
    -q Trif_boissieri \
    -q Trif_grandiflorum \
    -q Trif_strictum \
    -q Trif_aureum \
    -q Trif_glanduliferum \
    -q Trif_meduseum \
    -q Cic_arietinum \
    -o Cic_arietinum

# Vigna
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Vigna \
    -t Vig_angularis \
    -q Vig_unguiculata \
    -q Vig_radiata \
    -q Phas_vulgaris \
    -o Phas_vulgaris

#----------------------------#
# Fagaceae 壳斗科
#----------------------------#
# Quercus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Quercus \
    -t Q_rubra \
    -q Q_aliena \
    -q Q_spinosa \
    -q Q_aquifolioides \
    -q Trig_doichangensis \
    -o Trig_doichangensis

#----------------------------#
# Apocynaceae 夹竹桃科
#----------------------------#
# Asclepias
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Asclepias \
    -t Asc_syriaca \
    -q Asc_nivea \
    -q Cathar_roseus \
    -o Cathar_roseus

#----------------------------#
# Gentianaceae 龙胆科
#----------------------------#
# Erodium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Erodium \
    -t Ero_texanum \
    -q Ero_absinthoides \
    -q Ero_carvifolium \
    -q Ero_chrysanthum \
    -q Ero_crassifolium \
    -q Ero_gruinum \
    -q Ero_trifolium \
    -q Ger_palmatum \
    -o Ger_palmatum

# Pelargonium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pelargonium \
    -t Pel_x_hortorum \
    -q Pel_alternans \
    -q Ger_palmatum \
    -o Ger_palmatum

#----------------------------#
# Lentibulariaceae 狸藻科
#----------------------------#
# Utricularia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Utricularia \
    -t U_gibba \
    -q U_macrorhiza \
    -q Ping_ehlersiae \
    -o Ping_ehlersiae

#----------------------------#
# Oleaceae 木犀科
#----------------------------#
# Olea
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Olea \
    -t Olea_europaea \
    -q Olea_europaea_subsp_cuspidata \
    -q Olea_europaea_subsp_europaea \
    -q Olea_europaea_subsp_maroccana \
    -q Olea_woodiana_subsp_woodiana \
    -q Jas_nudiflorum \
    -o Jas_nudiflorum

#----------------------------#
# Orobanchaceae 列当科
#----------------------------#
# Cistanche
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cistanche \
    -t Cis_phelypaea \
    -q Cis_deserticola \
    -q Lin_philippensis \
    -o Lin_philippensis

# Orobanche
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Orobanche \
    -t Oro_ramosa \
    -q Oro_crenata \
    -q Oro_gracilis \
    -q Oro_purpurea \
    -q Oro_californica \
    -q Lin_philippensis \
    -o Lin_philippensis

# Fritillaria
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Fritillaria \
    -t Fri_cirrhosa \
    -q Fri_hupehensis \
    -q Fri_taipaiensis \
    -q Lil_superbum \
    -o Lil_superbum

#----------------------------#
# Melanthiaceae 黑药花科
#----------------------------#
# Trillium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Trillium \
    -t Tril_cuneatum \
    -q Tril_decumbens \
    -q Paris_verticillata \
    -o Paris_verticillata

#----------------------------#
# Magnoliaceae 木兰科
#----------------------------#
# Magnolia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Magnolia \
    -t Mag_grandiflora \
    -q Mag_salicifolia \
    -q Mag_pyramidata \
    -q Mag_tripetala \
    -q Mag_kobus \
    -q Mag_denudata \
    -q Mag_liliifera \
    -q Mag_officinalis \
    -q Mag_officinalis_subsp_biloba \
    -q Mag_kwangsiensis \
    -q Mag_cathcartii \
    -q Mag_dealbata \
    -q Mag_sprengeri \
    -q Mag_sinica \
    -q Mag_yunnanensis \
    -q Lir_tulipifera \
    -o Lir_tulipifera

#----------------------------#
# Chrysobalanaceae 金壳果科
#----------------------------#
# Hirtella
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Hirtella \
    -t Hir_physophora \
    -q Hir_racemosa \
    -q Chryso_icaco \
    -o Chryso_icaco

# Licania
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Licania \
    -t Lic_heteromorpha \
    -q Lic_alba \
    -q Lic_sprucei \
    -q Chryso_icaco \
    -o Chryso_icaco

#----------------------------#
# Salicaceae 杨柳科
#----------------------------#
# Populus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Populus \
    -t Pop_trichocarpa \
    -q Pop_alba \
    -q Pop_balsamifera \
    -q Pop_euphratica \
    -q Pop_tremula \
    -q Pop_yunnanensis \
    -q Pop_fremontii \
    -q Pop_cathayana \
    -q Salix_interior \
    -o Salix_interior

# Salix
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Salix \
    -t Salix_interior \
    -q Salix_purpurea \
    -q Salix_suchowensis \
    -q Pop_trichocarpa \
    -o Pop_trichocarpa

#----------------------------#
# Malvaceae 锦葵科
#----------------------------#
# Gossypium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Gossypium \
    -t Gos_barbadense \
    -q Gos_hirsutum \
    -q Gos_robinsonii \
    -q Gos_arboreum \
    -q Gos_raimondii \
    -q Gos_thurberi \
    -q Gos_herbaceum \
    -q Gos_mustelinum \
    -q Gos_darwinii \
    -q Gos_tomentosum \
    -q Gos_sturtianum \
    -q Gos_longicalyx \
    -q Gos_gossypioides \
    -q Gos_turneri \
    -q Gos_anomalum \
    -q Gos_stocksii \
    -q Gos_bickii \
    -q Gos_capitis_viridis \
    -q Gos_somalense \
    -q Gos_areysianum \
    -q Gos_incanum \
    -q Gos_herbaceum_subsp_africanum \
    -q The_cacao \
    -o The_cacao

#----------------------------#
# Myrtaceae 桃金娘科
#----------------------------#
# Angophora
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Angophora \
    -t Ango_floribunda \
    -q Ango_costata \
    -q Al_ternata \
    -o Al_ternata

# Corymbia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Corymbia \
    -t Cory_eximia \
    -q Cory_tessellaris \
    -q Cory_maculata \
    -q Cory_gummifera \
    -q Al_ternata \
    -o Al_ternata

# Eucalyptus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Eucalyptus \
    -t Eu_erythrocorys \
    -q Eu_microcorys \
    -q Eu_camaldulensis \
    -q Eu_curtisii \
    -q Eu_cloeziana \
    -q Eu_obliqua \
    -q Eu_grandis \
    -q Eu_globulus_subsp_globulus \
    -q Eu_aromaphloia \
    -q Eu_delegatensis \
    -q Eu_diversifolia \
    -q Eu_elata \
    -q Eu_nitens \
    -q Eu_radiata \
    -q Eu_regnans \
    -q Eu_sieberi \
    -q Eu_spathulata \
    -q Eu_umbra \
    -q Eu_saligna \
    -q Eu_deglupta \
    -q Eu_diversicolor \
    -q Eu_guilfoylei \
    -q Eu_marginata \
    -q Eu_melliodora \
    -q Eu_salmonophloia \
    -q Eu_torquata \
    -q Eu_cladocalyx \
    -q Eu_patens \
    -q Eu_verrucata \
    -q Eu_baxteri \
    -q Eu_polybractea \
    -q Al_ternata \
    -o Al_ternata

#----------------------------#
# Nymphaeaceae 睡莲科
#----------------------------#
# Nymphaea
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Nymphaea \
    -t Ny_alba \
    -q Ny_mexicana \
    -q Nu_advena \
    -o Nu_advena

#----------------------------#
# Poaceae 禾本科
#----------------------------#
# Aegilops
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Aegilops \
    -t Aeg_bicornis \
    -q Aeg_longissima \
    -q Aeg_searsii \
    -q Aeg_speltoides \
    -q Aeg_tauschii \
    -q Aeg_sharonensis \
    -q Aeg_cylindrica \
    -q Aeg_kotschyi \
    -q Aeg_geniculata \
    -q Zea_mays \
    -o Zea_mays

# Anthoxanthum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Anthoxanthum \
    -t Ant_odoratum \
    -q Ant_nitens \
    -q Zea_mays \
    -o Zea_mays

# Arundinaria
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Arundinaria \
    -t Aru_gigantea \
    -q Aru_tecta \
    -q Aru_fargesii \
    -q Aru_appalachiana \
    -q Zea_mays \
    -o Zea_mays

# Bambusa
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Bambusa \
    -t Bam_multiplex \
    -q Bam_oldhamii \
    -q Bam_bambos \
    -q Bam_emeiensis \
    -q Bam_arnhemica \
    -q Zea_mays \
    -o Zea_mays

# Chusquea
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Chusquea \
    -t Chu_circinata \
    -q Chu_liebmannii \
    -q Chu_spectabilis \
    -q Zea_mays \
    -o Zea_mays

# Fargesia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Fargesia \
    -t Far_yunnanensis \
    -q Far_nitida \
    -q Far_spathacea \
    -q Zea_mays \
    -o Zea_mays

# Festuca
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Festuca \
    -t Fes_arundinacea \
    -q Fes_pratensis \
    -q Fes_ovina \
    -q Fes_altissima \
    -q Zea_mays \
    -o Zea_mays

# Hordeum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Hordeum \
    -t Ho_jubatum \
    -q Ho_vulgare_subsp_vulgare \
    -q Zea_mays \
    -o Zea_mays

# Indocalamus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Indocalamus \
    -t Indoc_longiauritus \
    -q Indoc_wilsonii \
    -q Zea_mays \
    -o Zea_mays

# Lolium
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Lolium \
    -t Lol_multiflorum \
    -q Lol_perenne \
    -q Zea_mays \
    -o Zea_mays

# Melica
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Melica \
    -t Mel_mutica \
    -q Mel_subulata \
    -q Zea_mays \
    -o Zea_mays

# Oryza
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Oryza \
    -t Oryza_longistaminata \
    -q Oryza_rufipogon \
    -q Oryza_sativa \
    -q Oryza_australiensis \
    -q Oryza_officinalis \
    -q Oryza_nivara \
    -q Oryza_punctata \
    -q Oryza_glaberrima \
    -q Oryza_sativa_Indica_Group \
    -q Oryza_sativa_Japonica_Group \
    -q Oryza_glumipatula \
    -q Oryza_meridionalis \
    -q Oryza_barthii \
    -q Zea_mays \
    -o Zea_mays

# Pariana
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pariana \
    -t Paria_radiciflora \
    -q Paria_campestris \
    -q Zea_mays \
    -o Zea_mays

# Pharus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pharus \
    -t Phar_latifolius \
    -q Phar_lappulaceus \
    -q Zea_mays \
    -o Zea_mays

# Phyllostachys
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Phyllostachys \
    -t Phyl_edulis \
    -q Phyl_propinqua \
    -q Phyl_nigra_var_henonis \
    -q Phyl_sulphurea \
    -q Zea_mays \
    -o Zea_mays

# Saccharum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Saccharum \
    -t Sac_hybrid_cultivar_SP80_3280 \
    -q Sac_hybrid_cultivar_NCo_310 \
    -q Zea_mays \
    -o Zea_mays

# Sartidia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Sartidia \
    -t Sart_perrieri \
    -q Sart_dewinteri \
    -q Zea_mays \
    -o Zea_mays

# Sorghum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Sorghum \
    -t Sor_bicolor \
    -q Sor_timorense \
    -q Zea_mays \
    -o Zea_mays

# Triticum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Triticum \
    -t Trit_aestivum \
    -q Trit_monococcum \
    -q Trit_timopheevii \
    -q Trit_turgidum \
    -q Trit_urartu \
    -q Zea_mays \
    -o Zea_mays

#----------------------------#
# Rosaceae 蔷薇科
#----------------------------#
# Fragaria
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Fragaria \
    -t Fra_iinumae \
    -q Fra_chiloensis \
    -q Fra_virginiana \
    -q Fra_vesca_subsp_vesca \
    -q Fra_vesca_subsp_bracteata \
    -q Fra_mandshurica \
    -q Pri_utilis \
    -o Pri_utilis

# Prunus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Prunus \
    -t Pru_yedoensis \
    -q Pru_persica \
    -q Pru_maximowiczii \
    -q Pru_padus \
    -q Pru_mume \
    -q Pru_kansuensis \
    -q Pri_utilis \
    -o Pri_utilis

# Pyrus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pyrus \
    -t Py_pyrifolia \
    -q Py_spinosa \
    -q Pri_utilis \
    -o Pri_utilis

#----------------------------#
# Convolvulaceae 旋花科
#----------------------------#
# Cuscuta
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cuscuta \
    -t Cus_reflexa \
    -q Cus_gronovii \
    -q Cus_obtusiflora \
    -q Cus_exaltata \
    -q Ip_batatas \
    -o Ip_batatas

# Ipomoea
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Ipomoea \
    -t Ip_batatas \
    -q Ip_purpurea \
    -q Cus_reflexa \
    -o Cus_reflexa

#----------------------------#
# Solanaceae 茄科
#----------------------------#
# Capsicum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Capsicum \
    -t Capsi_annuum \
    -q Capsi_lycianthoides \
    -q Dat_stramonium \
    -o Dat_stramonium

# Dunalia
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Dunalia \
    -t Du_brachyacantha \
    -q Du_obovata \
    -q Du_solanacea \
    -q Dat_stramonium \
    -o Dat_stramonium

# Iochroma
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Iochroma \
    -t Io_loxense \
    -q Io_nitidum \
    -q Io_stenanthum \
    -q Io_tingoanum \
    -q Dat_stramonium \
    -o Dat_stramonium

# Nicotiana
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Nicotiana \
    -t Ni_sylvestris \
    -q Ni_tabacum \
    -q Ni_tomentosiformis \
    -q Ni_undulata \
    -q Dat_stramonium \
    -o Dat_stramonium

# Solanum
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Solanum \
    -t Sol_lycopersicum \
    -q Sol_peruvianum \
    -q Sol_chilense \
    -q Sol_pimpinellifolium \
    -q Sol_tuberosum \
    -q Sol_pennellii \
    -q Sol_habrochaites \
    -q Sol_neorickii \
    -q Sol_cheesmaniae \
    -q Sol_bulbocastanum \
    -q Sol_galapagense \
    -q Dat_stramonium \
    -o Dat_stramonium

#----------------------------------------------------------#
# Gymnosperm
#----------------------------------------------------------#
#----------------------------#
# Cupressaceae 柏科
#----------------------------#
# Callitropsis
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Callitropsis \
    -t Call_nootkatensis \
    -q Call_vietnamensis \
    -q Calo_formosana \
    -o Calo_formosana

# Juniperus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Juniperus \
    -t Ju_virginiana \
    -q Ju_scopulorum \
    -q Ju_monosperma \
    -q Ju_bermudiana \
    -q Calo_formosana \
    -o Calo_formosana

# Taiwania
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Taiwania \
    -t Tai_cryptomerioides \
    -q Tai_flousiana \
    -q Cun_lanceolata \
    -o Cun_lanceolata

#----------------------------#
# Taxaceae 红豆杉科 TODO
#----------------------------#
# Amentotaxus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Amentotaxus \
    -t Ame_argotaenia \
    -q Ame_formosana

# Cephalotaxus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cephalotaxus \
    -t Cep_wilsoniana \
    -q Cep_oliveri

#----------------------------#
# Podocarpaceae 罗汉松科
#----------------------------#
# Podocarpus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Podocarpus \
    -t Pod_totara \
    -q Pod_lambertii \
    -q Re_piresii \
    -o Re_piresii

#----------------------------#
# Cycadaceae
#----------------------------#
# Cycas
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Cycas \
    -t Cyc_revoluta \
    -q Cyc_taitungensis \
    -q Cyc_debaoensis \
    -q Bow_serrulata \
    -o Bow_serrulata

#----------------------------#
# Pinaceae 松科
#----------------------------#
# Picea
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Picea \
    -t Pic_abies \
    -q Pic_sitchensis \
    -q Pic_morrisonicola \
    -q Ced_deodara \
    -o Ced_deodara

# Pinus
perl ~/Scripts/withncbi/taxon/strain_bz.pl \
    --file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    --parallel 4 \
    --seq_dir ~/data/organelle/plastid_genomes \
    --use_name \
    --name Pinus \
    -t Pinus_contorta \
    -q Pinus_krempfii \
    -q Pinus_lambertiana \
    -q Pinus_strobus \
    -q Pinus_thunbergii \
    -q Pinus_taeda \
    -q Pinus_gerardiana \
    -q Pinus_nelsonii \
    -q Pinus_koraiensis \
    -q Pinus_massoniana \
    -q Pinus_monophylla \
    -q Pinus_taiwanensis \
    -q Ced_deodara \
    -o Ced_deodara
