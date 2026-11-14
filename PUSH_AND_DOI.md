# GitHub 푸시 → Zenodo DOI 발급 절차

이 저장소는 커밋까지 끝나 있습니다. 아래 순서대로 하시면 됩니다.

## 1. GitHub 저장소 생성 후 푸시

GitHub에서 **빈 저장소** `chelfit` 을 하나 만듭니다 (README·LICENSE 체크 해제 — 여기 이미 있습니다).

```bash
cd <이 폴더>
git remote add origin https://github.com/bisu9082/chelfit.git
git push -u origin main
```

## 2. Zenodo 연동 (최초 1회)

1. https://zenodo.org 에 GitHub 계정으로 로그인
2. 우상단 메뉴 → **GitHub**
3. 목록에서 `bisu9082/chelfit` 토글을 **ON**
   - 토글을 켠 *이후* 의 릴리스부터 DOI가 발급됩니다. 순서가 중요합니다.

## 3. 릴리스 발행 → DOI 자동 발급

```bash
git tag -a v1.0.0 -m "Version accompanying the submitted manuscript"
git push origin v1.0.0
```
그 다음 GitHub → Releases → **Draft a new release** → 태그 `v1.0.0` 선택 → Publish.

몇 분 안에 Zenodo가 DOI를 발급합니다.

## 4. 두 종류의 DOI — 원고에 넣을 것은 Concept DOI

Zenodo는 DOI를 두 개 줍니다.

| 종류 | 성격 | 원고 사용 |
|---|---|---|
| **Concept DOI** | 모든 버전을 대표. 항상 최신판으로 연결 | **이것을 쓰세요** |
| Version DOI | v1.0.0 한 판만 가리킴 | 특정 판을 못 박을 때만 |

Zenodo 레코드 페이지 오른쪽 **Versions** 칸에
*"Cite all versions? You can cite all versions by using the DOI ..."* 로 표시되는 것이 Concept DOI입니다.

## 5. 원고 3곳에 DOI 기입

`manuscript/` 안에서 아래 자리표시자를 바꾸세요.

| 파일 | 찾을 문자열 |
|---|---|
| `titleblock.tex` | `See DOI: \textit{[to be assigned]}` → ESI/저장소 DOI |
| `tail.tex` | `[repository DOI to be inserted at acceptance]` |
| `main.tex` | 위 둘이 합쳐진 파일이므로 **재조립 또는 직접 수정** |

재조립 명령:
```bash
cd manuscript
python3 - <<'PY'
pre=open('preamble.tex').read(); tb=open('titleblock.tex').read()
bd=open('body.tex').read(); fg=open('figures.tex').read(); tl=open('tail.tex').read()
bd=bd.replace('\\subsection{The effect under audit is large}', fg+'\n\\subsection{The effect under audit is large}',1)
open('main.tex','w').write(pre+tb+bd+tl)
PY
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## 6. 투고 전 마지막 확인

- [ ] **이름 영문 표기 확정** — 현재 `Jeongyun Kim`, `Jin Yoo`, `Ku Kang`, `Dongyeol Lee`.
      본인 확인 필수 (기존 논문 표기와 일치시킬 것)
- [ ] 이동열 선배 **이메일 주소** — `titleblock.tex` 의 `[to be supplied]`
- [ ] 육군사관학교 **정확한 영문 기관명·우편번호** (현재 `Korea Military Academy, Seoul 01805`)
- [ ] ORCID — `.zenodo.json` 의 각 creator 에 `"orcid": "0000-..."` 추가 가능
- [ ] Acknowledgements 작성
- [ ] Author contributions 각자 확인 (`tail.tex`)
