# 푸시 및 Zenodo DOI 발급 절차

이 저장소는 컨테이너 안에서 커밋만 되어 있습니다. 원격 푸시는 여기서 불가능하므로
아래 명령을 **본인 컴퓨터에서** 실행하십시오.

## 1. 저장소를 본인 컴퓨터로 가져오기

이 세션의 `chelfit-repo` 폴더를 통째로 내려받은 뒤, 그 폴더에서:

```bash
cd chelfit-repo
git log --oneline          # 커밋 3개 확인 (a126622가 최신)
```

## 2. GitHub에 빈 저장소 생성

https://github.com/new 에서 저장소를 만듭니다.
- 이름: `chelfit`
- 공개(Public)
- README / .gitignore / license **추가하지 않음** (이미 들어 있음)

## 3. 원격 연결 후 푸시

```bash
git remote add origin https://github.com/<본인계정>/chelfit.git
git branch -M main
git push -u origin main
```

`gh` CLI가 설치돼 있으면 2–3단계를 한 번에:

```bash
gh repo create chelfit --public --source=. --remote=origin --push
```

## 4. Zenodo 연동 및 DOI 발급

1. https://zenodo.org 에 GitHub 계정으로 로그인
2. 우측 상단 계정 → **GitHub**
3. 목록에서 `chelfit` 의 스위치를 **ON**
   (목록에 없으면 **Sync now** 를 먼저 누르십시오)
4. GitHub로 돌아와 릴리스 생성:

```bash
git tag -a v1.0 -m "Submission version for Chemical Science"
git push origin v1.0
gh release create v1.0 --title "v1.0 — submission" \
  --notes "Dataset, geometries, analysis code and manuscript accompanying the Chemical Science submission."
```

(`gh` 없이 하려면 GitHub 웹에서 Releases → Draft a new release → tag `v1.0` → Publish)

5. 릴리스가 게시되면 Zenodo가 자동으로 보관하고 DOI를 발급합니다.
   Zenodo의 **Upload** 목록에서 DOI를 확인하십시오. 두 개가 보입니다:
   - **Concept DOI** (모든 버전을 가리킴) ← 논문에 쓸 것
   - 버전 DOI (v1.0 전용)

## 5. 발급받은 DOI를 원고에 반영

`manuscript/tail.tex` 의 Data availability 절에서

```
[repository DOI to be inserted at acceptance]
```

를 실제 DOI로 바꾸고, `manuscript/titleblock.tex` 의 ESI 각주

```
See DOI: [to be assigned]
```

도 함께 교체한 뒤 다시 빌드하십시오:

```bash
cd manuscript
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## 투고 전 남은 항목

- [ ] 이동열 교수 이메일 주소 (`titleblock.tex` 의 `[to be supplied]`)
- [ ] Acknowledgements 본문 (`tail.tex`, 현재 `[To be completed.]`)
- [ ] Author contributions 각자 확인 (`tail.tex` 말미 대괄호 주석 삭제)
- [ ] Zenodo DOI 삽입 (위 5단계)
- [ ] ESI 파일(`si.tex`) 작성 — 현재 미작성
