# -*- coding: utf-8 -*-
"""원고의 모든 수치 ↔ JSON 산출물 전수 대조 [C-3]."""
import json,re,os
J={}
for f in ["experiment_summary.json","qm_summary_solv.json","perm_solv.json","rank_metric.json",
          "effect_summary.json","pair_summary.json","levels.json","physics_ladder.json",
          "pca_summary.json","power.json"]:
    if os.path.exists(f): J[f]=json.load(open(f))
def flat(d,p=""):
    o={}
    if isinstance(d,dict):
        for k,v in d.items(): o.update(flat(v,f"{p}.{k}"))
    elif isinstance(d,list):
        for i,v in enumerate(d): o.update(flat(v,f"{p}[{i}]"))
    elif isinstance(d,(int,float)): o[p]=float(d)
    return o
V={}
for f,d in J.items(): V.update(flat(d,f))
CLAIM=[
 ("중앙 배수",2.62,"effect"),("최대 배수",25.3,"effect"),("오차비율",0.110,"effect"),
 ("천장 R2",0.890,"effect"),("단일오차",0.0647,"effect"),("대비오차",0.0915,"effect"),
 ("LOPO 범용순위",-0.250,"model"),("LOPO 고전",-0.219,"model"),("LOPO QM용매",-0.340,"model"),
 ("LOPO 고전+QM",-0.480,"model"),
 ("순열 귀무평균",-0.207,"perm"),("순열 p(QM)",0.868,"perm"),("순열 p(고전)",0.766,"perm"),
 ("주효과 R2",0.133,"var"),("상호작용 R2",0.498,"var"),("in-sample QM용매",0.139,"var"),
 ("사다리 기체상복합",0.070,"ladder"),("사다리 기체상교환",0.086,"ladder"),
 ("사다리 용매복합",0.129,"ladder"),("사다리 용매교환",0.139,"ladder"),
 ("사다리 LOPO 기체복합",-0.227,"ladder"),("사다리 LOPO 용매교환",-0.340,"ladder"),
 ("풀링 r (용매)",0.227,"pair"),("풀링 p (용매)",0.027,"pair"),
 ("풀링 r (기체)",0.326,"pair"),("Spearman (기체)",0.241,"pair"),
 ("LOPO 크기(거리)",-0.277,"pair"),("LOPO 크기(벡터)",-0.240,"pair"),
 ("LOPO 방향 AUC",0.476,"pair"),("LOPO 부호",-0.563,"pair"),
 ("수준 킬레이터",-0.020,"level"),("수준 벡터",-0.088,"level"),("수준 논문",-0.086,"level"),
 ("수준 계열",-0.177,"level"),("수준 금속",-0.187,"level"),
 ("PC1",0.459,"pca"),("PC2",0.230,"pca"),("MW p",0.0008,"pca"),
]
bad=0
for name,c,_ in CLAIM:
    tol=max(6e-4,abs(c)*0.006)
    hit=[k for k,v in V.items() if abs(v-c)<tol]
    if hit: print(f"OK  {name:<20}{c:>9}  ← {hit[0]}")
    else:   print(f"!!! {name:<20}{c:>9}  ← 산출물에서 찾지 못함"); bad+=1
print(f"\n=== 대조 {len(CLAIM)}건 중 미확인 {bad}건 ===")
