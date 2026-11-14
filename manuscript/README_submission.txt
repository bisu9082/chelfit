chelfit — RSC Chemical Science Edge Article 투고 패키지
========================================================
빌드: pdflatex main → bibtex main → pdflatex main ×2

파일
  main.tex              원고 (RSC 2단 9pt, thpark_chemsci 템플릿 기반)
  main.pdf              조판 결과 (6면)
  cover_letter.tex/pdf  커버레터 (3면)
  references.bib        참고문헌 17건 — CrossRef 실데이터, RSC 양식
  rsc.bst               RSC BibTeX 스타일
  head_foot/            RSC 머리말·꼬리말 이미지
  Fig1–4.png            그림 4종 (400 dpi)
  preamble/titleblock/body/figures/tail .tex  조립용 분할 소스

템플릿 대비 변경
  - droidsans 미설치 환경 대비 조건부 로드로 교체
    \IfFileExists{droidsans.sty}{...}{\usepackage{helvet}}
    선생님 환경에 droidsans 가 있으면 자동으로 원래 폰트를 씁니다.
  - \wet / \wetnote 매크로 제거 (이 원고엔 불필요)
  - hyphenation 목록을 본 원고 용어로 교체
  - cover_letter.tex: letter 클래스에 없는 \paragraph 를 소제목용으로 정의

투고 전 확정 필요
  1. 공저자 명단 및 Author contributions (현재 Ku Kang 단독 + 각주 표시)
  2. Acknowledgements
  3. 데이터 저장소 DOI (Zenodo 등) → Data availability 와 ESI 각주에 기입
  4. ESI(si.tex) 작성 — 템플릿의 si.tex 구조 참고
