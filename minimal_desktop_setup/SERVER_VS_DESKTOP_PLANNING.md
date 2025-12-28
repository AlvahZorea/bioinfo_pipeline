# Server vs Desktop Testing - Planning Document

Use this document to plan what to run where, based on what you already have installed.

---

## Quick Reference: Resource Requirements

| Analysis | CPU | RAM | Database Size | Time (per sample) |
|----------|-----|-----|---------------|-------------------|
| **QC (FastQC/fastp)** | Low | 2 GB | None | 2-5 min |
| **Assembly (SPAdes)** | High | 16-32 GB | None | 30-60 min |
| **Assembly (Flye)** | Medium | 16 GB | None | 15-30 min |
| **Kraken2 (mini)** | Medium | 10 GB | 8 GB | 5-10 min |
| **Kraken2 (standard)** | Medium | 64 GB | 100 GB | 5-10 min |
| **Bakta (light)** | Medium | 8 GB | 2 GB | 10-20 min |
| **Bakta (full)** | Medium | 8 GB | 30 GB | 10-20 min |
| **GTDB-Tk** | High | 64 GB+ | 85 GB | 30-60 min |
| **AMRFinderPlus** | Low | 4 GB | 500 MB | 1-2 min |
| **MLST** | Low | 2 GB | 1 GB | <1 min |

---

## Inventory Checklist

### What's on your SERVER?

Fill in what you already have installed:

#### Tools Installed
- [ ] SPAdes version: _______
- [ ] Flye version: _______
- [ ] Unicycler version: _______
- [ ] Kraken2 version: _______
- [ ] GTDB-Tk version: _______
- [ ] Bakta version: _______
- [ ] Prokka version: _______
- [ ] AMRFinderPlus version: _______
- [ ] MLST version: _______
- [ ] ABRicate version: _______
- [ ] CheckV version: _______
- [ ] Pharokka version: _______
- [ ] Other: _______

#### Databases Available
- [ ] Kraken2 Standard (~100 GB): _______
- [ ] Kraken2 Mini (~8 GB): _______
- [ ] GTDB-Tk (~85 GB): _______
- [ ] Bakta Full (~30 GB): _______
- [ ] Bakta Light (~2 GB): _______
- [ ] CheckV (~2 GB): _______
- [ ] geNomad (~3 GB): _______
- [ ] CARD (~500 MB): _______
- [ ] AMRFinderPlus (~500 MB): _______
- [ ] VFDB (~100 MB): _______
- [ ] Pharokka (~8 GB): _______
- [ ] Other: _______

### What's on your DESKTOP/LAPTOP?

Specs:
- OS: _______
- RAM: _______ GB
- Free disk space: _______ GB
- CPU cores: _______

Tools to install:
- [ ] Conda/Mamba
- [ ] (list from minimal setup)

---

## Recommended Split: Server vs Desktop

Based on typical setups, here's what makes sense to run where:

### RUN ON SERVER (heavy computation, large DBs)

| Task | Why Server? |
|------|-------------|
| **Kraken2 with Standard DB** | Needs 64+ GB RAM |
| **GTDB-Tk classification** | Needs 64+ GB RAM, 85 GB database |
| **InterProScan** | 50 GB database |
| **Large assemblies** | Benefits from many cores |
| **Batch processing** | Multiple samples at once |

### RUN ON DESKTOP (lightweight, quick iteration)

| Task | Why Desktop? |
|------|-------------|
| **FastQC / MultiQC** | Quick, no database needed |
| **Small assemblies** | If you have 16+ GB RAM |
| **Bakta (light)** | Only needs 2 GB database |
| **AMRFinderPlus** | Small database, fast |
| **MLST** | Small database, instant |
| **Report generation** | Just data aggregation |
| **Pipeline development** | Quick iteration |

### EITHER (flexible)

| Task | Notes |
|------|-------|
| **Kraken2 with Mini DB** | 8 GB DB, needs ~10 GB RAM |
| **Assembly QC (QUAST)** | Low resources |
| **BLAST searches** | Depends on database size |

---

## Workflow Strategy Options

### Option A: Desktop Development, Server Validation
```
DESKTOP                          SERVER
   │                                │
   ├─ Develop pipeline              │
   ├─ Test with mini DBs            │
   ├─ Quick iterations              │
   │                                │
   └──── Transfer pipeline ────────►├─ Run with full DBs
                                    ├─ Validate results
                                    └─ Production runs
```

### Option B: Server Pre-computation, Desktop Assembly
```
SERVER                           DESKTOP
   │                                │
   ├─ Run Kraken2 (full)            │
   ├─ Run GTDB-Tk                   │
   ├─ Generate intermediate files   │
   │                                │
   └──── Transfer results ─────────►├─ Import pre-computed
                                    ├─ Run remaining analysis
                                    └─ Generate reports
```

### Option C: Full Server Pipeline (Desktop for viewing only)
```
SERVER                           DESKTOP
   │                                │
   ├─ Run entire pipeline           │
   ├─ Generate all outputs          │
   │                                │
   └──── Transfer reports ─────────►└─ View HTML reports
                                       Review results
```

---

## Test Samples Planning

### Sample 1: _______________________
- Type: [ ] Bacteria [ ] Virus [ ] Mixed
- Reads: [ ] Illumina [ ] Nanopore [ ] Both
- Expected species: _______
- Known features: _______
- Run where: [ ] Desktop [ ] Server [ ] Both

### Sample 2: _______________________
- Type: [ ] Bacteria [ ] Virus [ ] Mixed
- Reads: [ ] Illumina [ ] Nanopore [ ] Both
- Expected species: _______
- Known features: _______
- Run where: [ ] Desktop [ ] Server [ ] Both

### Sample 3: _______________________
- Type: [ ] Bacteria [ ] Virus [ ] Mixed
- Reads: [ ] Illumina [ ] Nanopore [ ] Both
- Expected species: _______
- Known features: _______
- Run where: [ ] Desktop [ ] Server [ ] Both

---

## Decision Matrix

Answer these questions to determine the best approach:

1. **How much RAM does your laptop have?**
   - < 8 GB → Server for most analysis
   - 8-16 GB → Desktop for light analysis
   - 16+ GB → Desktop can handle most things

2. **How much free disk space on laptop?**
   - < 20 GB → Minimal DBs only
   - 20-50 GB → Can add a few more DBs
   - 50+ GB → Comfortable for development

3. **Do you need offline/air-gapped capability?**
   - Yes → Everything must be local
   - No → Can use server for heavy tasks

4. **What's your iteration speed priority?**
   - Fast iteration → Desktop with mini DBs
   - Accuracy first → Server with full DBs

5. **What organisms are you testing?**
   - Common pathogens → Mini DBs probably fine
   - Rare/novel organisms → Need full DBs

---

## Notes

_Use this space for our discussion:_




