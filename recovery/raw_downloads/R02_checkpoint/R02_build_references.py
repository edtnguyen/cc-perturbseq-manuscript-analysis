from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path('/mnt/data')
SOURCE = ROOT / 'R00_source_snapshot.docx'
INVENTORY = ROOT / 'R01_citation_inventory.csv'
EXPECTED_SOURCE_SHA256 = '104f457aec0a2958871f13b4dbfa089d4f4fbbaea355e261665a3c7a1fb941ec'


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def rec(
    key: str,
    pmid: str,
    authors: list[str],
    title: str,
    journal: str,
    journal_abbrev: str,
    year: int,
    volume: str,
    issue: str,
    pages: str,
    doi: str,
    *,
    pmcid: str = '',
    publication_type: str = 'Research article',
) -> dict[str, Any]:
    return {
        'key': key,
        'pmid': pmid,
        'pmcid': pmcid,
        'authors': authors,
        'title': title,
        'journal': journal,
        'journal_abbrev': journal_abbrev,
        'year': str(year),
        'volume': volume,
        'issue': issue,
        'pages': pages,
        'doi': doi,
        'publication_type': publication_type,
        'pubmed_url': f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/',
        'doi_url': f'https://doi.org/{doi}',
    }


# Records are transcribed from the explicit PMIDs/author-year markers in the frozen manuscript
# and verified against PubMed and/or publisher metadata. Authors are retained in PubMed-style
# surname + initials form to avoid inventing expanded given names.
RECORDS: dict[str, dict[str, Any]] = {
    'Liberali2024Cell': rec(
        'Liberali2024Cell', '38906136', ['Liberali, P', 'Schier, AF'],
        'The evolution of developmental biology through conceptual and technological revolutions',
        'Cell', 'Cell', 2024, '187', '14', '3461-3495', '10.1016/j.cell.2024.05.053',
        publication_type='Review'),
    'Davidson2010Nature': rec(
        'Davidson2010Nature', '21164479', ['Davidson, EH'],
        'Emerging properties of animal gene regulatory networks',
        'Nature', 'Nature', 2010, '468', '7326', '911-920', '10.1038/nature09645',
        publication_type='Review'),
    'Ang2022Cell': rec(
        'Ang2022Cell', '35738284',
        ['Ang, LT', 'Nguyen, AT', 'Liu, KJ', 'Chen, A', 'Xiong, X', 'Curtis, M',
         'Martin, RM', 'Raftrey, BC', 'Ng, CY', 'Vogel, U', 'Lander, A', 'Lesch, BJ',
         'Fowler, JL', 'Holman, AR', 'Chai, T', 'Vijayakumar, S', 'Suchy, FP',
         'Nishimura, T', 'Bhadury, J', 'Porteus, MH', 'Nakauchi, H', 'Cheung, C',
         'George, SC', 'Red-Horse, K', 'Prescott, JB', 'Loh, KM'],
        'Generating human artery and vein cells from pluripotent stem cells highlights the arterial tropism of Nipah and Hendra viruses',
        'Cell', 'Cell', 2022, '185', '14', '2523-2541.e30', '10.1016/j.cell.2022.05.024',
        pmcid='PMC9617591'),
    'Shalaby1995Nature': rec(
        'Shalaby1995Nature', '7596435',
        ['Shalaby, F', 'Rossant, J', 'Gertsenstein, M', 'Yamaguchi, TP', 'Wu, XF', 'Breitman, ML', 'Schuh, AC'],
        'Failure of blood-island formation and vasculogenesis in Flk-1-deficient mice',
        'Nature', 'Nature', 1995, '376', '6535', '62-66', '10.1038/376062a0'),
    'Carmeliet1996Nature': rec(
        'Carmeliet1996Nature', '8602241',
        ['Carmeliet, P', 'Ferreira, V', 'Breier, G', 'Pollefeyt, S', 'Kieckens, L',
         'Gertsenstein, M', 'Fahrig, M', 'Vandenhoeck, A', 'Harpal, K', 'Eberhardt, C',
         'Declercq, C', 'Pawling, J', 'Moons, L', 'Collen, D', 'Risau, W', 'Nagy, A'],
        'Abnormal blood vessel development and lethality in embryos lacking a single VEGF allele',
        'Nature', 'Nature', 1996, '380', '6573', '435-439', '10.1038/380435a0'),
    'Ferrara1996Nature': rec(
        'Ferrara1996Nature', '8602242',
        ['Ferrara, N', 'Carver-Moore, K', 'Chen, H', 'Dowd, M', 'Lu, L', "O'Shea, KS", 'Powell-Braxton, L', 'Hillan, KJ', 'Moore, MW'],
        'Heterozygous embryonic lethality induced by targeted inactivation of the VEGF gene',
        'Nature', 'Nature', 1996, '380', '6573', '439-442', '10.1038/380439a0'),
    'Lawson2002DevCell': rec(
        'Lawson2002DevCell', '12110173', ['Lawson, ND', 'Vogel, AM', 'Weinstein, BM'],
        'Sonic hedgehog and vascular endothelial growth factor act upstream of the Notch pathway during arterial endothelial differentiation',
        'Developmental Cell', 'Dev. Cell', 2002, '3', '1', '127-136', '10.1016/S1534-5807(02)00198-3'),
    'Chetty2017DevBiol': rec(
        'Chetty2017DevBiol', '28279709',
        ['Chetty, SC', 'Rost, MS', 'Enriquez, JR', 'Schumacher, JA', 'Baltrunaite, K', 'Rossi, A', 'Stainier, DYR', 'Sumanas, S'],
        'Vegf signaling promotes vascular endothelial differentiation by modulating etv2 expression',
        'Developmental Biology', 'Dev. Biol.', 2017, '424', '2', '147-161', '10.1016/j.ydbio.2017.03.005'),
    'NussleinVolhard1980Nature': rec(
        'NussleinVolhard1980Nature', '6776413', ['Nüsslein-Volhard, C', 'Wieschaus, E'],
        'Mutations affecting segment number and polarity in Drosophila',
        'Nature', 'Nature', 1980, '287', '5785', '795-801', '10.1038/287795a0'),
    'Wieschaus2016ARCellDevBiol': rec(
        'Wieschaus2016ARCellDevBiol', '27501451', ['Wieschaus, E', 'Nüsslein-Volhard, C'],
        'The Heidelberg Screen for Pattern Mutants of Drosophila: A Personal Account',
        'Annual Review of Cell and Developmental Biology', 'Annu. Rev. Cell Dev. Biol.', 2016, '32', '', '1-46',
        '10.1146/annurev-cellbio-113015-023138', publication_type='Historical review'),
    'Brenner1974Genetics': rec(
        'Brenner1974Genetics', '4366476', ['Brenner, S'],
        'The genetics of Caenorhabditis elegans',
        'Genetics', 'Genetics', 1974, '77', '1', '71-94', '10.1093/genetics/77.1.71', pmcid='PMC1213120'),
    'Driever1996Development': rec(
        'Driever1996Development', '9007227',
        ['Driever, W', 'Solnica-Krezel, L', 'Schier, AF', 'Neuhauss, SC', 'Malicki, J', 'Stemple, DL',
         'Stainier, DY', 'Zwartkruis, F', 'Abdelilah, S', 'Rangini, Z', 'Belak, J', 'Boggs, C'],
        'A genetic screen for mutations affecting embryogenesis in zebrafish',
        'Development', 'Development', 1996, '123', '', '37-46', '10.1242/dev.123.1.37'),
    'Haffter1996Development': rec(
        'Haffter1996Development', '9007226',
        ['Haffter, P', 'Granato, M', 'Brand, M', 'Mullins, MC', 'Hammerschmidt, M', 'Kane, DA',
         'Odenthal, J', 'van Eeden, FJM', 'Jiang, YJ', 'Heisenberg, CP', 'Kelsh, RN',
         'Furutani-Seiki, M', 'Vogelsang, E', 'Beuchle, D', 'Schach, U', 'Fabian, C', 'Nüsslein-Volhard, C'],
        'The identification of genes with unique and essential functions in the development of the zebrafish, Danio rerio',
        'Development', 'Development', 1996, '123', '', '1-36', '10.1242/dev.123.1.1'),
    'Rosen2024NatCommun': rec(
        'Rosen2024NatCommun', '39419994',
        ['Rosen, BP', 'Li, QV', 'Cho, HS', 'Liu, D', 'Yang, D', 'Graff, S', 'Yan, J', 'Luo, R',
         'Verma, N', 'Damodaran, JR', 'Kale, HT', 'Kaplan, SJ', 'Beer, MA', 'Sidoli, S', 'Huangfu, D'],
        'Parallel genome-scale CRISPR-Cas9 screens uncouple human pluripotent stem cell identity versus fitness',
        'Nature Communications', 'Nat. Commun.', 2024, '15', '1', '8966', '10.1038/s41467-024-53284-4', pmcid='PMC11487130'),
    'Naxerova2021GenesDev': rec(
        'Naxerova2021GenesDev', '34711655',
        ['Naxerova, K', 'Di Stefano, B', 'Makofske, JL', 'Watson, EV', 'de Kort, MA', 'Martin, TD',
         'Dezfulian, M', 'Ricken, D', 'Wooten, EC', 'Kuroda, MI', 'Hochedlinger, K', 'Elledge, SJ'],
        'Integrated loss- and gain-of-function screens define a core network governing human embryonic stem cell behavior',
        'Genes & Development', 'Genes Dev.', 2021, '35', '21-22', '1527-1547', '10.1101/gad.349048.121', pmcid='PMC8559676'),
    'Li2019NatGenet': rec(
        'Li2019NatGenet', '31110351',
        ['Li, QV', 'Dixon, G', 'Verma, N', 'Rosen, BP', 'Gordillo, M', 'Luo, R', 'Xu, C', 'Wang, Q',
         'Soh, CL', 'Yang, D', 'Crespo, M', 'Shukla, A', 'Xiang, Q', 'Dündar, F', 'Zumbo, P',
         'Witkin, M', 'Koche, R', 'Betel, D', 'Chen, S', 'Massagué, J', 'Garippa, R', 'Evans, T',
         'Beer, MA', 'Huangfu, D'],
        'Genome-scale screens identify JNK-JUN signaling as a barrier for pluripotency exit and endoderm differentiation',
        'Nature Genetics', 'Nat. Genet.', 2019, '51', '6', '999-1010', '10.1038/s41588-019-0408-9', pmcid='PMC6545159'),
    'Xu2025NatCellBiol': rec(
        'Xu2025NatCellBiol', '40897805',
        ['Xu, N', 'Cho, HS', 'Hackland, JOS', 'Benito-Kwiecinski, S', 'Saurat, N', 'Harschnitz, O',
         'Russo, MV', 'Garippa, R', 'Ciceri, G', 'Studer, L'],
        'Genome-wide CRISPR screen identifies Menin and SUZ12 as regulators of human developmental timing',
        'Nature Cell Biology', 'Nat. Cell Biol.', 2025, '27', '9', '1411-1421', '10.1038/s41556-025-01751-5'),
    'Xu2020StemCells': rec(
        'Xu2020StemCells', '32129551',
        ['Xu, J', 'Zhou, C', 'Foo, KS', 'Yang, R', 'Xiao, Y', 'Bylund, K', 'Sahara, M', 'Chien, KR'],
        'Genome-wide CRISPR screen identifies ZIC2 as an essential gene that controls the cell fate of early mesodermal precursors to human heart progenitors',
        'Stem Cells', 'Stem Cells', 2020, '38', '6', '741-755', '10.1002/stem.3168'),
    'Padmanabhan2024NatCardiovascRes': rec(
        'Padmanabhan2024NatCardiovascRes', '39196112',
        ['Padmanabhan, A', 'de Soysa, TY', 'Pelonero, A', 'Sapp, V', 'Shah, PP', 'Wang, Q', 'Li, L',
         'Lee, CY', 'Sadagopan, N', 'Nishino, T', 'Ye, L', 'Yang, R', 'Karnay, A', 'Poleshko, A',
         'Bolar, N', 'Linares-Saldana, R', 'Ranade, SS', 'Alexanian, M', 'Morton, SU', 'Jain, M',
         'Haldar, SM', 'Srivastava, D', 'Jain, R'],
        'A genome-wide CRISPR screen identifies BRD4 as a regulator of cardiomyocyte differentiation',
        'Nature Cardiovascular Research', 'Nat. Cardiovasc. Res.', 2024, '3', '3', '317-331', '10.1038/s44161-024-00431-1'),
    'Lee2024Circulation': rec(
        'Lee2024Circulation', '38752370',
        ['Lee, CJM', 'Autio, MI', 'Zheng, W', 'Song, Y', 'Wang, SC', 'Wong, DCP', 'Xiao, J',
         'Zhu, Y', 'Yusoff, P', 'Chock, WK', 'Low, BC', 'Sudol, M', 'Foo, RSY'],
        'Genome-Wide CRISPR Screen Identifies an NF2-Adherens Junction Mechanistic Dependency for Cardiac Lineage',
        'Circulation', 'Circulation', 2024, '149', '25', '1960-1979', '10.1161/CIRCULATIONAHA.122.061335'),
    'Dixit2016Cell': rec(
        'Dixit2016Cell', '27984732',
        ['Dixit, A', 'Parnas, O', 'Li, B', 'Chen, J', 'Fulco, CP', 'Jerby-Arnon, L', 'Marjanovic, ND',
         'Dionne, D', 'Burks, T', 'Raychowdhury, R', 'Adamson, B', 'Norman, TM', 'Lander, ES',
         'Weissman, JS', 'Friedman, N', 'Regev, A'],
        'Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens',
        'Cell', 'Cell', 2016, '167', '7', '1853-1866.e17', '10.1016/j.cell.2016.11.038'),
    'Adamson2016Cell': rec(
        'Adamson2016Cell', '27984733',
        ['Adamson, B', 'Norman, TM', 'Jost, M', 'Cho, MY', 'Núñez, JK', 'Chen, Y', 'Villalta, JE',
         'Gilbert, LA', 'Horlbeck, MA', 'Hein, MY', 'Pak, RA', 'Gray, AN', 'Gross, CA', 'Dixit, A',
         'Parnas, O', 'Regev, A', 'Weissman, JS'],
        'A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response',
        'Cell', 'Cell', 2016, '167', '7', '1867-1882.e21', '10.1016/j.cell.2016.11.048', pmcid='PMC5315571'),
    'Jaitin2016Cell': rec(
        'Jaitin2016Cell', '27984734',
        ['Jaitin, DA', 'Weiner, A', 'Yofe, I', 'Lara-Astiaso, D', 'Keren-Shaul, H', 'David, E',
         'Salame, TM', 'Tanay, A', 'van Oudenaarden, A', 'Amit, I'],
        'Dissecting Immune Circuits by Linking CRISPR-Pooled Screens with Single-Cell RNA-Seq',
        'Cell', 'Cell', 2016, '167', '7', '1883-1896.e15', '10.1016/j.cell.2016.11.039'),
    'Kotliar2019eLife': rec(
        'Kotliar2019eLife', '31282856',
        ['Kotliar, D', 'Veres, A', 'Nagy, MA', 'Tabrizi, S', 'Hodis, E', 'Melton, DA', 'Sabeti, PC'],
        'Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq',
        'eLife', 'eLife', 2019, '8', '', 'e43803', '10.7554/eLife.43803'),
    'Schnitzler2024Nature': rec(
        'Schnitzler2024Nature', '38326615',
        ['Schnitzler, GR', 'Kang, H', 'Fang, S', 'Angom, RS', 'Lee-Kim, VS', 'Ma, XR', 'Zhou, R',
         'Zeng, T', 'Guo, K', 'Taylor, MS', 'Vellarikkal, SK', 'Barry, AE', 'Sias-Garcia, O',
         'Bloemendal, A', 'Munson, G', 'Guckelberger, P', 'Nguyen, TH', 'Bergman, DT',
         'Hinshaw, S', 'Cheng, N', 'Cleary, B', 'Aragam, K', 'Lander, ES', 'Finucane, HK',
         'Mukhopadhyay, D', 'Gupta, RM', 'Engreitz, JM'],
        'Convergence of coronary artery disease genes onto endothelial cell programs',
        'Nature', 'Nature', 2024, '626', '8000', '799-807', '10.1038/s41586-024-07022-x', pmcid='PMC10921916'),
    'Takahashi2001EMBOJ': rec(
        'Takahashi2001EMBOJ', '11387210', ['Takahashi, T', 'Yamaguchi, S', 'Chida, K', 'Shibuya, M'],
        'A single autophosphorylation site on KDR/Flk-1 is essential for VEGF-A-dependent activation of PLC-gamma and DNA synthesis in vascular endothelial cells',
        'The EMBO Journal', 'EMBO J.', 2001, '20', '11', '2768-2778', '10.1093/emboj/20.11.2768', pmcid='PMC125481'),
    'Nishizuka1992Science': rec(
        'Nishizuka1992Science', '1411571', ['Nishizuka, Y'],
        'Intracellular signaling by hydrolysis of phospholipids and activation of protein kinase C',
        'Science', 'Science', 1992, '258', '5082', '607-614', '10.1126/science.1411571', publication_type='Review'),
    'Ma2024AdvSci': rec(
        'Ma2024AdvSci', '38161226',
        ['Ma, T', 'Liu, X', 'Su, H', 'Shi, Q', 'He, Y', 'Wu, F', 'Gao, C', 'Li, K', 'Liang, Z',
         'Zhang, D', 'Zhang, X', 'Hu, K', 'Li, S', 'Wang, L', 'Wang, M', 'Yue, S', 'Hong, W',
         'Chen, X', 'Zhang, J', 'Zheng, L', 'Deng, X', 'Wang, P', 'Fan, Y'],
        'Coupling of Perinuclear Actin Cap and Nuclear Mechanics in Regulating Flow-Induced Yap Spatiotemporal Nucleocytoplasmic Transport',
        'Advanced Science', 'Adv. Sci. (Weinh.)', 2024, '11', '11', 'e2305867', '10.1002/advs.202305867', pmcid='PMC10953556'),
    'Cowdin2025ATVB': rec(
        'Cowdin2025ATVB', '40501385',
        ['Cowdin, MA', 'Pramanik, T', 'Mohr-Allen, SR', 'Fu, Y', 'Mills, A', 'Spurgin, SB',
         'Varner, VD', 'Davis, GE', 'Cleaver, O'],
        'Lats1/2 Are Essential for Developmental Vascular Remodeling and Biomechanical Adaptation to Shear Stress',
        'Arteriosclerosis, Thrombosis, and Vascular Biology', 'Arterioscler. Thromb. Vasc. Biol.',
        2025, '45', '9', '1521-1542', '10.1161/ATVBAHA.124.322258', pmcid='PMC12648833'),
    'Nakajima2017DevCell': rec(
        'Nakajima2017DevCell', '28350986',
        ['Nakajima, H', 'Yamamoto, K', 'Agarwala, S', 'Terai, K', 'Fukui, H', 'Fukuhara, S',
         'Ando, K', 'Miyazaki, T', 'Yokota, Y', 'Schmelzer, E', 'Belting, HG', 'Affolter, M',
         'Lecaudey, V', 'Mochizuki, N'],
        'Flow-Dependent Endothelial YAP Regulation Contributes to Vessel Maintenance',
        'Developmental Cell', 'Dev. Cell', 2017, '40', '6', '523-536.e6', '10.1016/j.devcel.2017.02.019'),
    'Lundin2020DevCell': rec(
        'Lundin2020DevCell', '32032546',
        ['Lundin, V', 'Sugden, WW', 'Theodore, LN', 'Sousa, PM', 'Han, A', 'Chou, S', 'Wrighton, PJ',
         'Cox, AG', 'Ingber, DE', 'Goessling, W', 'Daley, GQ', 'North, TE'],
        'YAP Regulates Hematopoietic Stem Cell Formation in Response to the Biomechanical Forces of Blood Flow',
        'Developmental Cell', 'Dev. Cell', 2020, '52', '4', '446-460.e5', '10.1016/j.devcel.2020.01.006', pmcid='PMC7398148'),
    'Wang2016PNAS': rec(
        'Wang2016PNAS', '27671657',
        ['Wang, KC', 'Yeh, YT', 'Nguyen, P', 'Limqueco, E', 'Lopez, J', 'Thorossian, S',
         'Guan, KL', 'Li, YSJ', 'Chien, S'],
        'Flow-dependent YAP/TAZ activities regulate endothelial phenotypes and atherosclerosis',
        'Proceedings of the National Academy of Sciences of the United States of America',
        'Proc. Natl Acad. Sci. USA', 2016, '113', '41', '11525-11530', '10.1073/pnas.1613121113', pmcid='PMC5068257'),
    'Nimmo2013DevCell': rec(
        'Nimmo2013DevCell', '23911199',
        ['Nimmo, R', 'Ciau-Uitz, A', 'Ruiz-Herguido, C', 'Soneji, S', 'Bigas, A', 'Patient, R', 'Enver, T'],
        'miR-142-3p controls the specification of definitive hemangioblasts during ontogeny',
        'Developmental Cell', 'Dev. Cell', 2013, '26', '3', '237-249', '10.1016/j.devcel.2013.06.023'),
    'Loh2024SeminCellDevBiol': rec(
        'Loh2024SeminCellDevBiol', '37393122', ['Loh, KM', 'Ang, LT'],
        'Building human artery and vein endothelial cells from pluripotent stem cells, and enduring mysteries surrounding arteriovenous development',
        'Seminars in Cell & Developmental Biology', 'Semin. Cell Dev. Biol.', 2024, '155', 'Pt C', '62-75', '10.1016/j.semcdb.2023.06.004',
        publication_type='Review'),
}

# Explicit occurrence-to-record mapping from the R01 ledger.
CITATION_MAP = {
    'CIT-0001': 'Liberali2024Cell',
    'CIT-0002': 'Davidson2010Nature',
    'CIT-0003': 'Ang2022Cell',
    'CIT-0004': 'Shalaby1995Nature',
    'CIT-0005': 'Carmeliet1996Nature',
    'CIT-0006': 'Ferrara1996Nature',
    'CIT-0007': 'Lawson2002DevCell',
    'CIT-0008': 'Chetty2017DevBiol',
    'CIT-0009': 'Ang2022Cell',
    'CIT-0011': 'NussleinVolhard1980Nature',
    'CIT-0012': 'Wieschaus2016ARCellDevBiol',
    'CIT-0013': 'Brenner1974Genetics',
    'CIT-0014': 'Driever1996Development',
    'CIT-0015': 'Haffter1996Development',
    'CIT-0016': 'Rosen2024NatCommun',
    'CIT-0017': 'Naxerova2021GenesDev',
    'CIT-0018': 'Li2019NatGenet',
    'CIT-0019': 'Xu2025NatCellBiol',
    'CIT-0020': 'Xu2020StemCells',
    'CIT-0021': 'Padmanabhan2024NatCardiovascRes',
    'CIT-0022': 'Lee2024Circulation',
    'CIT-0023': 'Dixit2016Cell',
    'CIT-0024': 'Adamson2016Cell',
    'CIT-0025': 'Jaitin2016Cell',
    'CIT-0027': 'Kotliar2019eLife',
    'CIT-0030': 'Ang2022Cell',
    'CIT-0035': 'Schnitzler2024Nature',
    'CIT-0043': 'Takahashi2001EMBOJ',
    'CIT-0044': 'Nishizuka1992Science',
    'CIT-0046': 'Ma2024AdvSci',
    'CIT-0047': 'Cowdin2025ATVB',
    'CIT-0048': 'Nakajima2017DevCell',
    'CIT-0049': 'Lundin2020DevCell',
    'CIT-0050': 'Wang2016PNAS',
    'CIT-0057': 'Nimmo2013DevCell',
    'CIT-0058': 'Loh2024SeminCellDevBiol',
}

BATCH_B_IDS = {'CIT-0046', 'CIT-0047', 'CIT-0048', 'CIT-0049', 'CIT-0050', 'CIT-0057', 'CIT-0058'}

SUPPORT = {
    'CIT-0001': ('Contextual', 'Broad review supports the developmental-biology framing; it is not evidence for a specific regulatory map.'),
    'CIT-0002': ('Contextual', 'Review directly supports the concept of gene-regulatory networks coordinating development.'),
    'CIT-0003': ('Direct', 'Primary study directly describes rapid hPSC-to-arterial/venous endothelial differentiation through the stated intermediates.'),
    'CIT-0004': ('Foundational', 'Primary loss-of-function study establishes FLK1/KDR requirement for blood-island formation and vasculogenesis; it is not a protocol citation by itself.'),
    'CIT-0005': ('Foundational', 'Primary Vegf haploinsufficiency study establishes VEGF requirement for embryonic vascular development.'),
    'CIT-0006': ('Foundational', 'Independent primary Vegf haploinsufficiency study establishes VEGF requirement for embryonic vascular development.'),
    'CIT-0007': ('Direct', 'Primary zebrafish study directly links VEGF signaling to arterial endothelial differentiation.'),
    'CIT-0008': ('Direct', 'Primary study directly links Vegf signaling to vascular endothelial differentiation through etv2 expression.'),
    'CIT-0009': ('Duplicate exact match', 'Author-year marker is the same Cell 2022 record as PMID 35738284.'),
    'CIT-0011': ('Direct historical', 'Foundational forward genetic screen in Drosophila.'),
    'CIT-0012': ('Contextual historical', 'Retrospective account of the Heidelberg screen; appropriate as historical context, not independent primary evidence.'),
    'CIT-0013': ('Direct historical', 'Foundational forward genetic screen in C. elegans.'),
    'CIT-0014': ('Direct historical', 'Primary zebrafish embryogenesis mutagenesis screen.'),
    'CIT-0015': ('Direct historical', 'Primary zebrafish mutagenesis screen identifying essential developmental genes.'),
    'CIT-0016': ('Direct example', 'Genome-scale CRISPR screens in human pluripotent stem cells; supports the human-cell screening claim.'),
    'CIT-0017': ('Direct example', 'Loss- and gain-of-function screens in human embryonic stem cells.'),
    'CIT-0018': ('Direct example', 'Genome-scale screens during pluripotency exit and endoderm differentiation.'),
    'CIT-0019': ('Direct example', 'Genome-wide CRISPR screen of human developmental timing.'),
    'CIT-0020': ('Direct example', 'Genome-wide CRISPR screen in early mesoderm-to-cardiac differentiation.'),
    'CIT-0021': ('Direct example', 'Genome-wide CRISPR screen of cardiomyocyte differentiation.'),
    'CIT-0022': ('Direct example', 'Genome-wide CRISPR screen identifying an NF2-dependent cardiac-lineage mechanism.'),
    'CIT-0023': ('Direct method', 'Foundational Perturb-seq paper.'),
    'CIT-0024': ('Direct method', 'Foundational pooled single-cell CRISPR screening paper.'),
    'CIT-0025': ('Direct method', 'Foundational pooled CRISPR plus single-cell RNA-seq paper.'),
    'CIT-0027': ('Direct method', 'Original cNMF method paper.'),
    'CIT-0030': ('Duplicate exact match', 'Author-year marker is the same Cell 2022 record as PMID 35738284.'),
    'CIT-0035': ('Direct methodological precedent', 'Exact author-year match; prior endothelial Perturb-seq study used cNMF to define gene programs. Exact downstream testing details should still be described from the current study.'),
    'CIT-0043': ('Direct', 'Primary study directly establishes KDR/Flk-1-dependent PLC-gamma activation.'),
    'CIT-0044': ('Foundational', 'Canonical review supports phospholipid hydrolysis, diacylglycerol and PKC activation; it is not endothelial-specific.'),
    'CIT-0046': ('Direct', 'Primary study directly examines flow-shear-induced YAP transport in endothelial cells.'),
    'CIT-0047': ('Direct', 'Primary study directly tests LATS1/2-dependent endothelial adaptation to shear stress and vascular remodeling.'),
    'CIT-0048': ('Direct', 'Primary study directly links flow-dependent endothelial YAP regulation to vessel maintenance.'),
    'CIT-0049': ('Partial/contextual', 'Primary study concerns blood-flow-dependent YAP in hemogenic endothelium and hematopoietic stem-cell formation, not mature endothelial homeostasis.'),
    'CIT-0050': ('Direct', 'Primary study directly links flow-dependent endothelial YAP/TAZ activity to endothelial phenotypes and atherosclerosis.'),
    'CIT-0057': ('Partial', 'Primary study addresses definitive hemangioblast specification but does not by itself establish the complete HAND1→ETV2→TAL1/SCL temporal sequence stated in the manuscript.'),
    'CIT-0058': ('Contextual review', 'Review is relevant to human artery/vein endothelial differentiation, but a primary reference is preferable for the exact sequential-expression claim.'),
}


def split_pages(pages: str) -> tuple[str, str]:
    # Keep article numbers or electronic locators in SP only.
    if re.fullmatch(r'e?\d+[A-Za-z0-9.]*', pages):
        return pages, ''
    if '-' in pages:
        a, b = pages.split('-', 1)
        return a, b
    return pages, ''


def write_ris(path: Path, records: list[dict[str, Any]]) -> None:
    lines: list[str] = []
    for r in records:
        lines.append('TY  - JOUR')
        lines.append(f'ID  - {r["key"]}')
        for author in r['authors']:
            lines.append(f'AU  - {author}')
        lines.append(f'TI  - {r["title"]}')
        lines.append(f'JF  - {r["journal"]}')
        lines.append(f'JO  - {r["journal_abbrev"]}')
        lines.append(f'PY  - {r["year"]}')
        if r['volume']:
            lines.append(f'VL  - {r["volume"]}')
        if r['issue']:
            lines.append(f'IS  - {r["issue"]}')
        sp, ep = split_pages(r['pages'])
        if sp:
            lines.append(f'SP  - {sp}')
        if ep:
            lines.append(f'EP  - {ep}')
        lines.append(f'DO  - {r["doi"]}')
        lines.append(f'UR  - {r["doi_url"]}')
        lines.append(f'UR  - {r["pubmed_url"]}')
        lines.append(f'N1  - PMID: {r["pmid"]}')
        if r['pmcid']:
            lines.append(f'N1  - PMCID: {r["pmcid"]}')
        lines.append('ER  -')
        lines.append('')
    path.write_text('\n'.join(lines), encoding='utf-8')


def bib_escape(text: str) -> str:
    return text.replace('\\', '\\\\').replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')


def write_bib(path: Path, records: list[dict[str, Any]]) -> None:
    chunks: list[str] = []
    for r in records:
        fields = [
            ('author', ' and '.join(r['authors'])),
            ('title', '{' + r['title'] + '}'),
            ('journal', r['journal']),
            ('year', r['year']),
            ('volume', r['volume']),
            ('number', r['issue']),
            ('pages', r['pages'].replace('-', '--')),
            ('doi', r['doi']),
            ('pmid', r['pmid']),
            ('pmcid', r['pmcid']),
            ('url', r['doi_url']),
        ]
        lines = [f'@article{{{r["key"]},']
        for name, value in fields:
            if value:
                lines.append(f'  {name} = {{{bib_escape(value)}}},')
        lines[-1] = lines[-1].rstrip(',')
        lines.append('}')
        chunks.append('\n'.join(lines))
    path.write_text('\n\n'.join(chunks) + '\n', encoding='utf-8')


def nature_preview(r: dict[str, Any], number: int) -> str:
    authors = r['authors']
    if len(authors) <= 5:
        au = ', '.join(a.replace(', ', ' ') for a in authors[:-1])
        if len(authors) > 1:
            au += ' & ' + authors[-1].replace(', ', ' ')
        elif authors:
            au = authors[0].replace(', ', ' ')
    else:
        au = authors[0].replace(', ', ' ') + ' et al.'
    issue = f'({r["issue"]})' if r['issue'] else ''
    return f'{number}. {au} {r["title"]}. {r["journal_abbrev"]} {r["volume"]}{issue}, {r["pages"]} ({r["year"]}). https://doi.org/{r["doi"]}'


def main() -> None:
    assert SOURCE.exists() and SOURCE.stat().st_size > 0
    assert INVENTORY.exists() and INVENTORY.stat().st_size > 0
    assert sha256(SOURCE) == EXPECTED_SOURCE_SHA256, 'R00 source hash mismatch'

    with INVENTORY.open(newline='', encoding='utf-8-sig') as f:
        all_rows = list(csv.DictReader(f))
    r02_rows = [row for row in all_rows if row['planned_unit'] == 'R02']
    assert len(r02_rows) == 36
    ids = [row['citation_id'] for row in r02_rows]
    assert len(ids) == len(set(ids))
    assert set(ids) == set(CITATION_MAP)
    assert set(CITATION_MAP.values()) <= set(RECORDS)
    assert len(RECORDS) == 34

    # First-appearance order is used only as a provisional preview/map; final numbering occurs in R05.
    ordered_keys: list[str] = []
    for row in r02_rows:
        key = CITATION_MAP[row['citation_id']]
        if key not in ordered_keys:
            ordered_keys.append(key)
    assert len(ordered_keys) == 34

    a_rows = [r for r in r02_rows if r['citation_id'] not in BATCH_B_IDS]
    b_rows = [r for r in r02_rows if r['citation_id'] in BATCH_B_IDS]
    assert len(a_rows) == 29 and len(b_rows) == 7
    a_keys = [k for k in ordered_keys if any(CITATION_MAP[x['citation_id']] == k for x in a_rows)]
    b_keys = [k for k in ordered_keys if any(CITATION_MAP[x['citation_id']] == k for x in b_rows)]
    assert len(a_keys) == 27 and len(b_keys) == 7

    # Output paths.
    out = {
        'input_subset': ROOT / 'R02_input_inventory.csv',
        'a_ris': ROOT / 'R02A_existing_refs_intro_results.ris',
        'a_bib': ROOT / 'R02A_existing_refs_intro_results.bib',
        'a_audit': ROOT / 'R02A_metadata_audit.csv',
        'b_ris': ROOT / 'R02B_existing_refs_hippo_discussion.ris',
        'b_bib': ROOT / 'R02B_existing_refs_hippo_discussion.bib',
        'b_audit': ROOT / 'R02B_metadata_audit.csv',
        'all_ris': ROOT / 'R02_existing_references_all.ris',
        'all_bib': ROOT / 'R02_existing_references_all.bib',
        'map': ROOT / 'R02_reference_map_provisional.csv',
        'flags': ROOT / 'R02_support_flags.csv',
        'unresolved': ROOT / 'R02_unresolved_references.csv',
        'preview': ROOT / 'R02_Nature_reference_preview.txt',
        'report': ROOT / 'R02_resolution_report.txt',
        'manifest': ROOT / 'R02_checkpoint_manifest.txt',
        'validation': ROOT / 'R02_validation.json',
        'checksums': ROOT / 'R02_checksums.txt',
        'zip': ROOT / 'R02_checkpoint.zip',
        'zip_sha': ROOT / 'R02_checkpoint_sha256.txt',
    }

    # Copy the exact R02 input slice.
    with out['input_subset'].open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(r02_rows[0].keys()))
        writer.writeheader()
        writer.writerows(r02_rows)

    write_ris(out['a_ris'], [RECORDS[k] for k in a_keys])
    write_bib(out['a_bib'], [RECORDS[k] for k in a_keys])
    write_ris(out['b_ris'], [RECORDS[k] for k in b_keys])
    write_bib(out['b_bib'], [RECORDS[k] for k in b_keys])
    write_ris(out['all_ris'], [RECORDS[k] for k in ordered_keys])
    write_bib(out['all_bib'], [RECORDS[k] for k in ordered_keys])

    audit_fields = [
        'batch', 'citation_id', 'manuscript_page', 'paragraph_locator', 'section', 'subsection',
        'current_marker', 'record_key', 'resolution_status', 'resolution_basis', 'duplicate_of_record_key',
        'title', 'authors', 'journal', 'journal_abbrev', 'year', 'volume', 'issue', 'pages_or_article',
        'doi', 'pmid', 'pmcid', 'publication_type', 'bibliographic_confidence', 'support_level',
        'support_assessment', 'pubmed_url', 'doi_url', 'manuscript_claim_context', 'notes'
    ]

    def audit_row(row: dict[str, str], batch: str) -> dict[str, str]:
        key = CITATION_MAP[row['citation_id']]
        r = RECORDS[key]
        current = row['current_marker']
        if current.startswith('PMID:'):
            basis = 'Exact PMID supplied in manuscript'
            duplicate = ''
        elif current.startswith('Ang et al.'):
            basis = 'Exact author-year match to PMID 35738284; deduplicated with the PMID occurrence'
            duplicate = 'Ang2022Cell'
        elif current.startswith('Schnitzler et al.'):
            basis = 'Exact author-year/title/context match to PMID 38326615 and DOI 10.1038/s41586-024-07022-x'
            duplicate = ''
        else:
            basis = 'Resolved from manuscript marker and context'
            duplicate = ''
        support_level, assessment = SUPPORT[row['citation_id']]
        notes = ''
        if row['citation_id'] in {'CIT-0009', 'CIT-0030'}:
            notes = 'Do not create an additional bibliography entry; both author-year markers map to the same Cell 2022 record.'
        if row['citation_id'] == 'CIT-0048':
            notes = 'PMID 28350986 resolves to Nakajima et al., not the separate Kim et al. JCI angiogenesis paper (PMID 28805663).'
        return {
            'batch': batch,
            'citation_id': row['citation_id'],
            'manuscript_page': row['page'],
            'paragraph_locator': row['paragraph_locator'],
            'section': row['section_level_1'],
            'subsection': row['section_level_2'],
            'current_marker': current,
            'record_key': key,
            'resolution_status': 'Resolved',
            'resolution_basis': basis,
            'duplicate_of_record_key': duplicate,
            'title': r['title'],
            'authors': '; '.join(r['authors']),
            'journal': r['journal'],
            'journal_abbrev': r['journal_abbrev'],
            'year': r['year'],
            'volume': r['volume'],
            'issue': r['issue'],
            'pages_or_article': r['pages'],
            'doi': r['doi'],
            'pmid': r['pmid'],
            'pmcid': r['pmcid'],
            'publication_type': r['publication_type'],
            'bibliographic_confidence': 'High',
            'support_level': support_level,
            'support_assessment': assessment,
            'pubmed_url': r['pubmed_url'],
            'doi_url': r['doi_url'],
            'manuscript_claim_context': row['exact_claim_or_context'],
            'notes': notes,
        }

    a_audit_rows = [audit_row(r, 'R02A') for r in a_rows]
    b_audit_rows = [audit_row(r, 'R02B') for r in b_rows]
    for path, rows in [(out['a_audit'], a_audit_rows), (out['b_audit'], b_audit_rows)]:
        with path.open('w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=audit_fields)
            writer.writeheader()
            writer.writerows(rows)

    map_fields = ['provisional_number', 'record_key', 'first_citation_id', 'first_page', 'title', 'first_author', 'year', 'doi', 'pmid', 'occurrence_count', 'all_citation_ids']
    map_rows = []
    for n, key in enumerate(ordered_keys, 1):
        occurrences = [r for r in r02_rows if CITATION_MAP[r['citation_id']] == key]
        r = RECORDS[key]
        map_rows.append({
            'provisional_number': n,
            'record_key': key,
            'first_citation_id': occurrences[0]['citation_id'],
            'first_page': occurrences[0]['page'],
            'title': r['title'],
            'first_author': r['authors'][0],
            'year': r['year'],
            'doi': r['doi'],
            'pmid': r['pmid'],
            'occurrence_count': len(occurrences),
            'all_citation_ids': '; '.join(x['citation_id'] for x in occurrences),
        })
    with out['map'].open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=map_fields)
        writer.writeheader()
        writer.writerows(map_rows)

    flags = [
        {
            'flag_id': 'R02-F01', 'citation_ids': 'CIT-0046; CIT-0047; CIT-0048; CIT-0049; CIT-0050',
            'page': '7', 'paragraph_locator': 'P0084', 'severity': 'Needs wording or new-reference review',
            'issue': 'The five listed records support blood flow/shear-stress regulation of endothelial Hippo-YAP/TAZ signaling. None was verified as a direct hydrostatic-pressure study.',
            'recommended_R02_disposition': 'Keep all five records resolved, but do not treat them as sufficient support for the phrase “and hydrostatic pressure”. Queue a pressure-specific search for R07C or narrow the wording to flow/shear stress.'
        },
        {
            'flag_id': 'R02-F02', 'citation_ids': 'CIT-0049', 'page': '7', 'paragraph_locator': 'P0084',
            'severity': 'Contextual only',
            'issue': 'PMID 32032546 studies YAP in blood-flow-dependent hematopoietic stem-cell formation from hemogenic endothelium, not mature endothelial homeostasis.',
            'recommended_R02_disposition': 'Retain as developmental/biomechanical context only; do not present it as a mature-EC shear-response study.'
        },
        {
            'flag_id': 'R02-F03', 'citation_ids': 'CIT-0057', 'page': '9', 'paragraph_locator': 'P0112',
            'severity': 'Claim-reference mismatch',
            'issue': 'PMID 23911199 concerns miR-142-3p control of definitive hemangioblast specification and does not by itself establish the full HAND1→ETV2→TAL1/SCL sequence.',
            'recommended_R02_disposition': 'Retain only as related developmental context; search for direct primary evidence for the sequence in R07F or revise the sentence.'
        },
        {
            'flag_id': 'R02-F04', 'citation_ids': 'CIT-0058', 'page': '9', 'paragraph_locator': 'P0112',
            'severity': 'Review/contextual',
            'issue': 'PMID 37393122 is a review. It is relevant but should not be the sole support for an exact temporal sequence if a primary paper is available.',
            'recommended_R02_disposition': 'Retain the review for synthesis and add direct primary evidence during R07F.'
        },
        {
            'flag_id': 'R02-F05', 'citation_ids': 'CIT-0048', 'page': '7', 'paragraph_locator': 'P0084',
            'severity': 'Identifier disambiguation resolved',
            'issue': 'The manuscript PMID 28350986 resolves to Nakajima et al., “Flow-Dependent Endothelial YAP Regulation Contributes to Vessel Maintenance” (Dev. Cell 2017), not Kim et al., “YAP/TAZ regulates sprouting angiogenesis and vascular barrier maturation” (JCI 2017; PMID 28805663).',
            'recommended_R02_disposition': 'Use the Nakajima record for PMID 28350986. Consider the Kim paper only later if a separate claim needs it.'
        },
    ]
    with out['flags'].open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(flags[0].keys()))
        writer.writeheader()
        writer.writerows(flags)

    unresolved_fields = ['citation_id', 'current_marker', 'page', 'reason_unresolved', 'required_author_input', 'next_action']
    with out['unresolved'].open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=unresolved_fields)
        writer.writeheader()
        # Deliberately zero unresolved bibliographic identities in R02.

    out['preview'].write_text(
        'R02 PROVISIONAL NATURE-STYLE REFERENCE PREVIEW\n'
        'This is a metadata-checking preview only. Final numbering and Word insertion occur in R05.\n\n' +
        '\n'.join(nature_preview(RECORDS[key], i) for i, key in enumerate(ordered_keys, 1)) + '\n',
        encoding='utf-8')

    counts = Counter(CITATION_MAP.values())
    report = f'''R02 REFERENCE-RESOLUTION REPORT\n\nSTATUS: PASS\n\nSCOPE\n- Canonical source: R00_source_snapshot.docx\n- Canonical source SHA-256: {sha256(SOURCE)}\n- Input ledger: R01_citation_inventory.csv\n- Scope restricted to references already present in the Abstract, Introduction, Results and active Discussion.\n- No manuscript text, citation markers, comments, figures or formatting were modified.\n- No new literature was added for previously uncited claims in this unit.\n\nRESULTS\n- R02 occurrence-level markers reviewed: {len(r02_rows)}\n- Unique bibliographic records resolved: {len(ordered_keys)}\n- Exact PMID occurrences: {sum(1 for r in r02_rows if r['current_marker'].startswith('PMID:'))}\n- Author-year occurrences: {sum(1 for r in r02_rows if not r['current_marker'].startswith('PMID:'))}\n- Duplicate author-year occurrences merged into an existing PMID record: 2\n- Bibliographic identities unresolved: 0\n- R02A occurrences / unique records: {len(a_rows)} / {len(a_keys)}\n- R02B occurrences / unique records: {len(b_rows)} / {len(b_keys)}\n- Existing explicit citations found in active Discussion: 0\n\nIMPORTANT SUPPORT FLAGS\n1. The five Hippo mechanosignaling papers support endothelial flow/shear biology, but the listed set did not verify the manuscript's separate phrase “hydrostatic pressure”.\n2. PMID 32032546 is developmental hemogenic-endothelium/HSC context, not a mature-EC homeostasis paper.\n3. PMID 23911199 does not by itself establish the complete HAND1→ETV2→TAL1/SCL sequence.\n4. PMID 37393122 is a review; primary support should be added for the exact sequence if retained.\n5. PMID 28350986 was disambiguated as Nakajima et al., Dev. Cell 2017 (DOI 10.1016/j.devcel.2017.02.019). It is not the separate Kim et al. JCI paper (PMID 28805663).\n\nOUTPUTS\n- R02A_existing_refs_intro_results.ris / .bib: {len(a_keys)} unique records\n- R02A_metadata_audit.csv: {len(a_rows)} occurrence rows\n- R02B_existing_refs_hippo_discussion.ris / .bib: {len(b_keys)} unique records\n- R02B_metadata_audit.csv: {len(b_rows)} occurrence rows\n- R02_existing_references_all.ris / .bib: {len(ordered_keys)} unique merged records\n- R02_reference_map_provisional.csv: first-appearance map for later R05 insertion\n- R02_support_flags.csv: claim-support cautions\n- R02_unresolved_references.csv: header only; zero unresolved identities\n- R02_Nature_reference_preview.txt: provisional display-format check\n\nNEXT UNIT\nR03 resolves the references already present in Methods and Supplementary Information.\n'''
    out['report'].write_text(report, encoding='utf-8')

    # Validation before packaging.
    def count_ris(path: Path) -> int:
        return sum(1 for line in path.read_text(encoding='utf-8').splitlines() if line == 'TY  - JOUR')

    def count_bib(path: Path) -> int:
        return len(re.findall(r'^@article\{', path.read_text(encoding='utf-8'), flags=re.M))

    validation = {
        'status': 'PASS',
        'source': {
            'path': str(SOURCE),
            'sha256_expected': EXPECTED_SOURCE_SHA256,
            'sha256_observed': sha256(SOURCE),
            'unchanged': sha256(SOURCE) == EXPECTED_SOURCE_SHA256,
        },
        'inventory': {
            'r02_occurrences': len(r02_rows),
            'unique_citation_ids': len(set(ids)),
            'mapped_citation_ids': len(CITATION_MAP),
            'all_ids_mapped_exactly_once': set(ids) == set(CITATION_MAP),
        },
        'records': {
            'unique_records': len(ordered_keys),
            'unique_pmids': len({RECORDS[k]['pmid'] for k in ordered_keys}),
            'unique_dois_casefolded': len({RECORDS[k]['doi'].casefold() for k in ordered_keys}),
            'r02a_unique': len(a_keys),
            'r02b_unique': len(b_keys),
            'unresolved_identities': 0,
        },
        'formats': {
            'r02a_ris_entries': count_ris(out['a_ris']),
            'r02a_bib_entries': count_bib(out['a_bib']),
            'r02b_ris_entries': count_ris(out['b_ris']),
            'r02b_bib_entries': count_bib(out['b_bib']),
            'merged_ris_entries': count_ris(out['all_ris']),
            'merged_bib_entries': count_bib(out['all_bib']),
        },
        'support_flags': len(flags),
        'checks': [],
    }
    checks = [
        (validation['source']['unchanged'], 'R00 source hash unchanged'),
        (validation['inventory']['all_ids_mapped_exactly_once'], 'All 36 R02 citation IDs mapped exactly once'),
        (len(ordered_keys) == 34, '34 unique records after deduplication'),
        (validation['records']['unique_pmids'] == 34, 'No duplicate PMID records'),
        (validation['records']['unique_dois_casefolded'] == 34, 'No duplicate DOI records'),
        (count_ris(out['a_ris']) == count_bib(out['a_bib']) == 27, 'R02A RIS/BibTeX parity'),
        (count_ris(out['b_ris']) == count_bib(out['b_bib']) == 7, 'R02B RIS/BibTeX parity'),
        (count_ris(out['all_ris']) == count_bib(out['all_bib']) == 34, 'Merged RIS/BibTeX parity'),
        (len(a_audit_rows) == 29 and len(b_audit_rows) == 7, 'Occurrence audit row counts'),
        (all(RECORDS[k]['title'] and RECORDS[k]['authors'] and RECORDS[k]['journal'] and RECORDS[k]['year'] and RECORDS[k]['doi'] and RECORDS[k]['pmid'] for k in ordered_keys), 'Required metadata present'),
    ]
    validation['checks'] = [{'name': name, 'passed': bool(passed)} for passed, name in checks]
    if not all(passed for passed, _ in checks):
        validation['status'] = 'FAIL'
    out['validation'].write_text(json.dumps(validation, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    assert validation['status'] == 'PASS', validation

    # Manifest and checksums for all files except package and package checksum.
    member_keys = [
        'input_subset', 'a_ris', 'a_bib', 'a_audit', 'b_ris', 'b_bib', 'b_audit',
        'all_ris', 'all_bib', 'map', 'flags', 'unresolved', 'preview', 'report',
        'validation'
    ]
    manifest = f'''R02 CHECKPOINT MANIFEST\n\nSTATUS: PASS\nCanonical source: R00_source_snapshot.docx\nCanonical source SHA-256: {sha256(SOURCE)}\nR02 occurrence rows: {len(r02_rows)}\nUnique resolved records: {len(ordered_keys)}\nR02A: {len(a_rows)} occurrences / {len(a_keys)} records\nR02B: {len(b_rows)} occurrences / {len(b_keys)} records\nBibliographic identities unresolved: 0\nSupport flags retained for later review: {len(flags)}\nManuscript modified in R02: NO\n\nFILES\n'''
    for key in member_keys:
        p = out[key]
        manifest += f'- {p.name}: {p.stat().st_size:,} bytes; SHA-256 {sha256(p)}\n'
    manifest += f'- {Path(__file__).name}: {Path(__file__).stat().st_size:,} bytes; SHA-256 {sha256(Path(__file__))}\n'
    out['manifest'].write_text(manifest, encoding='utf-8')

    checksum_members = [out[k] for k in member_keys] + [out['manifest'], Path(__file__)]
    out['checksums'].write_text(
        ''.join(f'{sha256(p)}  {p.name}\n' for p in checksum_members), encoding='utf-8')

    package_members = checksum_members + [out['checksums']]
    if out['zip'].exists():
        out['zip'].unlink()
    with zipfile.ZipFile(out['zip'], 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for p in package_members:
            z.write(p, arcname=p.name)
    with zipfile.ZipFile(out['zip']) as z:
        assert z.testzip() is None
        assert set(z.namelist()) == {p.name for p in package_members}
    out['zip_sha'].write_text(f'{sha256(out["zip"])}  {out["zip"].name}\n', encoding='utf-8')

    # Final existence and nonzero checks.
    final_paths = [out[k] for k in out]
    for p in final_paths:
        assert p.exists() and p.stat().st_size > 0, p

    print(json.dumps({
        'status': 'PASS',
        'source_sha256': sha256(SOURCE),
        'r02_occurrences': len(r02_rows),
        'unique_records': len(ordered_keys),
        'r02a_occurrences': len(a_rows),
        'r02a_records': len(a_keys),
        'r02b_occurrences': len(b_rows),
        'r02b_records': len(b_keys),
        'support_flags': len(flags),
        'checkpoint_zip': str(out['zip']),
        'checkpoint_zip_bytes': out['zip'].stat().st_size,
        'checkpoint_zip_sha256': sha256(out['zip']),
    }, indent=2))


if __name__ == '__main__':
    main()
