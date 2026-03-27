# Quick NR Simulation Test Run Lines: 

## Python scripts: 

Parallel processing logging -> chart generation
```python
python3 -m log_to_excel_v2 -d logs/ -ex . -sn parallel_500_slots

python3 -m make_charts -e runlog_results.xlsx -c charts.yml -o img
```

## Checking output progress: 

```bash
du -sh logs/

ls logs/ | wc -l 
```



# Quick Git Reference: 

### Standard git flow from feat/ branch
```git
git status
git add . 
git commit -m ""
git push -u origin feat/branch_name
```

### Merge branch -> Main 
```git
git status
git checkout main 
git pull origin main 
git status
git merge feat/branch_name
```

### Merge Main -> branch
```git
git checkout branch_name
git status
git fetch origin 
git merge origin/main 
```