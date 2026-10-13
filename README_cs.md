# 🌐 Language / Jazyk
* [English](README.md)
* [Čeština](#čeština)

---

# Čeština
# Diferenciace druhů rodu *Penicillium* (série *Viridicata*) za použití MALDI-TOF MS a strojového učení

Tento repozitář obsahuje bioinformatický skript vyvinutý během zahraniční výzkumné stáže Erasmus+ ve **Westerdijk Fungal Biodiversity Institute** (Utrecht, Nizozemsko). Projekt je zaměřen na automatizovanou druhovou diferenciaci uvnitř taxonomicky náročné série *Penicillium* *Viridicata* s využitím spektrálních dat z MALDI-TOF hmotnostní spektrometrie a pokročilých statistických modelů v jazyce R.

---

## Přehled projektu

### Problém
Tradiční morfologická identifikace blízce příbuzných druhů rodu *Penicillium* je extrémně časově náročná a vyžaduje vysokou expertízu. Sekvenování DNA sice poskytuje přesné výsledky, ale pro rutinní diagnostiku s vysokým objemem vzorků je drahé a pomalé.

### Řešení
Hmotnostní spektrometrie MALDI-TOF nabízí rychlou, přesnou a nákladově efektivní alternativu založenou na jedinečných proteinových profilech ("otiscích prstů"). Komerční softwary (jako uzavřený MBT Compass Explorer od Brukeru) však mají fixně nastavené parametry a často selhávají při čistém rozlišení vizuálně podobných druhových komplexů, což vede k nepřehledným a roztříštěným dendrogramům.

### Tento projekt
K překonání limitů komerčního softwaru byl vytvořen tento vlastní bioinformatický kód v R za využití open-source knihoven. Ořezáním surových spekter na optimální sektory, implementací striktnějšího zarovnání píků a nasazením algoritmů strojového učení se podařilo drasticky snížit chybovost klasifikace (OOB error rate) ze **40,68 %** na pouhých **8,7 %**.

---

## Metodika a pracovní postup

1. **Mikrobiologie a příprava vzorků:** 
   * Kultivace kmenů ze sbírky CBS/DTO na MEA agaru a médiu Sabouraud Broth.
   * Vysokokapacitní proteinová extrakce za využití kyseliny mravenčí (HCOOH) a acetonitrilu.
2. **Akvizice spekter a kontrola kvality:**
   * Automatické měření v programu `flexControl` (Bruker) – generování 32 spekter na jeden vzorek.
   * Manuální filtrace plochých spekter a odlehlých hodnot ve `flexAnalysis` (dosažení ≥ 20 kvalitních spekter na kmen).
3. **Bioinformatické zpracování (Skript v R):**
   * **Preprocessing:** Ořezání spektrálního rozsahu (2500–4700 Da), normalizace pomocí celkového iontového proudu (TIC) a průměrování spekter.
   * **Extrakce znaků:** Detekce píků (metoda MAD, SNR = 3) a striktní binování píků (tolerance = 0,0009).
   * **Statistická analýza:** Kruskal-Wallisův test s Benjamini-Hochbergovou korekcí falešné pozitivity pro identifikaci druhově specifických biomarkerů.
   * **Strojové učení a vícerozměrná statistika:** Hierarchické shlukování (Ward.D2), PCA, t-SNE a klasifikace pomocí Random Forest (500 stromů).

---

## Ochrana dat a anonymizace (Důležité)

Z důvodu ochrany duševního vlastnictví a interních pravidel Westerdijk institutu **nejsou v tomto repozitáři dostupná surová (raw) spektrální data**. 

Aby byla pipeline plně funkční a spustitelná pro demonstrační účely, je v repozitáři přiložen **matematicky upravený a plně anonymizovaný ukázkový dataset** (`data_matrix.rds`), ze kterého byla vypreparována veškerá skrytá metadata:
* Obsahuje oříznutý vzorek dat (8 modelových skupin označených anonymně jako `sample1` až `sample8`).
* Všechny identifikátory byly z objektu smazány.
* Samotné hodnoty intenzit píků byly matematicky upraveny (zašuměny), aby nebylo možné profily zpětně identifikovat oproti referenčním databázím, avšak geometrie distribuce dat pro výpočty v R zůstala plně zachována.

Vzhledem k absenci raw spekter je úvodní část skriptu `penicillium_maldi_analysis.R` (import, normalizace a detekce píků) zakomentována. Skript začíná přímo od načtení vyčištěné matice (je však třeba načíst R balíčky a nastavit pracovní adresář):

```R
data_matrix <- readRDS("data_matrix.rds")
```

---

## Struktura repozitáře

* `R/` - Složka obsahující zdrojový kód a zpracovaná data.
  * `data_matrix.rds` - Předzpracovaná datová matice připravená pro analýzu.
  * `penicillium_maldi_analysis.R` - Kompletní R skript obsahující preprocessing dat, identifikaci biomarkerů, PCA, t-SNE a Random Forest modelování.
* `LICENSE` - Licenční podmínky pro použití kódu.
* `README.md` - Dokumentace projektu v angličtině.
* `README_cs.md` - Dokumentace projektu v češtině.

---

## Systémové požadavky a R balíčky

### Verze R
* **R (v4.4.3)** – Tato specifická verze byla použita k zajištění plné zpětné kompatibility a stability s potřebnými bioinformatickými a hmotnostně-spektrometrickými balíčky.

### Potřebné knihovny
Pro instalaci všech balíčků nutných pro zpracování dat, strojové učení a vizualizaci spusťte v R konzoli následující příkaz:

```R
install.packages(c("readBrukerFlexData", "MALDIquant", "MALDIquantForeign", 
                   "MALDIrppa", "pheatmap", "Rtsne", "randomForest", 
                   "ggplot2", "ggrepel"))
```