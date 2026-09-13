# Dog Y-haplogroup nomenclature: Oetjens/Ding HG names vs Embark letters

Stage 9d calls haplogroups on the Oetjens et al. 2018 tree (BMC Genomics
19:350), which uses the Ding et al. 2012 (Heredity 108:507) names. Embark's
letter scheme (A1a, A1b, A2a, A2b, B, B1, C, D, E, F) is Embark's own
relabeling, introduced in a 2017 SMBE poster (Lounsberry, Khan, Wells, Fulop,
Sams, Boyko), layered on the Shannon et al. 2015 (PNAS) haplotype names.
Oetjens assigned haplogroups "using the definitions from Shannon et al." and
then split HG27 out of HG1-3, so both schemes share a haplotype vocabulary
(H1a, H3, H5a, H7, H8, H9, Ha, Hb, Hc, H27 ...) and map at that level.

| Oetjens HG | Embark | Confidence | Evidence |
|---|---|---|---|
| HG1-3 | A = A1a + A1b + A2a + A2b | certain for A1a; likely-to-certain for all of A | Poster Fig. 1: A1a = H1a.x, A1b = Ha.x, A2b = Hc.1/Hc.9-17 + H3, A2a = Hc.2-8. Ding: H1 (= H1a) in HG1, H3 in HG3; Oetjens lumps HG1+HG3 (77/118 canids, every Western breed). Anchor: Cosmo = Embark A1a/H1a.59 = our HG1-3. One-to-many: A1a/A1b/A2a/A2b need Shannon/Embark array markers our SNV tree lacks. |
| HG27 | C | certain (identical haplotype set) | Oetjens defines HG27 = Hb.1, H27, H5b, H5a; poster C bar = exactly those. Breeds: Jindo, Shiba, Akita, Tibetan Mastiff; Embark C blurb: Akita, Shiba, NGSD, Samoyed, Malamute. |
| HG6 | B (and very likely B1) | likely | Poster B = H15.x; Ding places H15 in HG6 (East Asia only). B1 (arctic, basal, wolf-introgressed per Embark) is the sister branch; whether B1 carries the HG6 diagnostic alleles is not demonstrated. Oetjens HG6: Tibetan Mastiffs, Chinese/Indian village dogs, one Yellowstone wolf. |
| HG23 | D | likely-to-certain | Poster D = H7, H10, Hd; Ding puts H7 (Afghan) and H10 (Saluki) in HG23. Oetjens HG23: Afghan x3, Saluki, Sloughi, Tibetan Terrier x2, village dogs, two wolves. Embark D: Afghan, Lhasa Apso, Sarplaninac, French Bulldog, Bull Terrier. |
| HG8 | E | certain | Poster E = H8.x; Oetjens HG8 = one Nigerian village dog, sister to HG23; Ding H8 = Basenji, Canaan Dog, Turkey, Thai Ridgeback. |
| HG9 | F | certain | Poster F = H9; Ding HG9 = H9 (Basenjis) + H18 (E. Siberian Laika); Oetjens HG9 = Nigerian village dog + Xinjiang wolf; Embark F "closer to wolves than to other dogs". |

Embark 2017 sample frequencies (mostly Western breeds): A1a 57.8, A2b 13.1,
A1b 8.5, A2a 2.2 (A total 81.6), D 11.9, C 3.2, B 2.1, B1 0.4, F 0.5, E 0.3 %.

Haplotype names: Ding H1-H31 (indel-aware H1a, H5a/b ...); Shannon 2015
subdivided with a dotted index (H1a.8) and coined Ha/Hb/Hc for unmatched
haplotypes; Embark added new indices (H1a.59, Hc.17, B1a-c, Hd.1). So Embark
"H1a.59" = Ding H1/H1a sub-haplotype 59, inside Embark A1a = Ding HG1 =
Oetjens HG1-3. Embark "C1" is a MITOCHONDRIAL label, not paternal.

Sources: Oetjens 2018 PMC5946424 (+ Additional file 1 sample table); Ding
2012 PMC3330686; Shannon 2015 PMC4640804 (Fig. 3); Brown 2011 PMC3237445;
Sacks 2013 MBE 30:1103; Embark/Boyko SMBE 2017 poster
http://front.embarkvet.com/conference-presentations/y_chromosome_variation_05_2017.pdf;
Embark blog "Embarking on Dog Ancestry Research" (Sams 2018); sled-dog paper
PMC8454093 Suppl. Table 1 (Embark letter -> haplotype lists).

App side: dogs-app/lib/yHaplogroups.ts carries the map + customer-facing lore.
