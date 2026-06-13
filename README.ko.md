# eHT-PCR

[English](README.md) | **한국어**

**한 좌위(locus)**만 증폭하고, 다른 위치는 증폭하지 않는 PCR 프라이머를 설계합니다.

eHT-PCR은 전산상에서 대량 PCR을 수행하기 위한 명령줄 도구입니다. **참조 서열**과
**타겟**을 입력하면 후보 프라이머 쌍을 생성하고, 각 프라이머를 참조 서열에 다시 정렬한 뒤,
PCR 산물이 **하나만** 나오는 쌍만 남깁니다. 즉 paralog, homeolog, 그 밖의 off-target을
함께 증폭할 가능성이 있는 프라이머 쌍을 걸러냅니다. 특정 생물종에 묶여 있지 않으며,
참조 서열은 유전체, 전사체, 또는 임의의 multi-FASTA일 수 있습니다.

## 동작 방식

각 타겟에 대해 eHT-PCR은 다음 순서로 동작합니다.

1. 타겟 서열을 **추출**합니다.
2. 길이, GC 함량, Tm, GC-clamp, 산물 크기 조건을 만족하는 후보 프라이머 쌍을 **생성**합니다.
3. `bwa`로 모든 프라이머를 특이성 검사 참조 서열에 **매핑**한 뒤, 정방향/역방향 hit를 짝지어 증폭 산물(amplicon)을 계산합니다.
4. on-target 증폭 산물이 하나만 나오는 프라이머 쌍만 **보존**합니다.

여기서 "unique"는 증폭 산물이 *하나*라는 뜻입니다. 프라이머 하나가 개별 조건을 잘
만족한다는 뜻이 아니라, 쌍으로 PCR을 했을 때 산물이 하나만 예측된다는 뜻입니다.

## 설치

제공되는 conda 환경은 eHT-PCR과 함께 사용할 수 있는 외부 도구(`bwa`, BLAST+, miniprot)를
같이 설치합니다.

```bash
conda env create -f environment.yml
conda activate ehtpcr
ehtpcr --help
```

이미 준비된 환경이 있다면 패키지만 설치할 수도 있습니다. Python 3.10 이상이 필요합니다.

```bash
pip install .
```

Tm 모델은 별도 설정 없이 바로 사용할 수 있습니다. nearest-neighbor 모델과 `primer3`는
Python 의존성으로 함께 설치됩니다. 다만 정렬 단계는 `PATH`에 있는 외부 실행 파일을
호출합니다. 특이성 검사에는 **`bwa`**가 필요합니다. `--no-spec`로 특이성 검사를 끄고
프라이머 설계만 수행하는 경우에는 필요하지 않습니다. 서열로 타겟 위치를 찾는
`--query-fasta` 모드에서는 **BLAST+** 및/또는 **miniprot**이 추가로 필요할 수 있습니다.

```bash
conda install -c bioconda bwa blast miniprot
```

위 명령으로 세 도구를 모두 설치할 수 있습니다.

## 빠른 시작: 전체 실행 예제

eHT-PCR에는 두 유전자가 들어 있는 작은 예제 참조 서열 `examples/mini.fa`가 포함되어
있습니다.

```text
>geneX example gene geneX
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGTCCGCAACACT
CGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGAAGACAACCCTCTGGTTCTTT
...
>geneY example gene geneY
CTTCTTTAGGCGAGAGTACCCTATTTTTGGCCCTATGAGCGCCTTGATGGACTCGTTACTTGGGACCAAT
...
```

함께 제공되는 `config/example.yaml`은 FASTA header에 `gene`이 들어 있는 모든 항목, 즉
두 유전자 모두에 대해 프라이머를 설계합니다. 각 프라이머 쌍의 특이성은 같은 파일을 기준으로
검사합니다.

```bash
ehtpcr run --config config/example.yaml --out result/
```

`result/candidates.tsv`에는 후보 프라이머 쌍이 한 행씩 기록됩니다. 실제로 사용할 후보는
**`kept`** 행입니다. 이 행들은 타겟에서 단일 증폭 산물(`n_valid_amplicons == 1`)을 내는
것으로 판정된 프라이머 쌍입니다. 아래는 일부 열만 표시한 예입니다.

```text
target_id  pair_id             f_seq                    r_seq                  product_size  n_valid_amplicons  unique  keep_decision
geneX      geneX:3-26:132-153  CTGCTGTCGGACTCCTAGTTACG  TCTTACGGACGGGAAAGAACC  150           1                  True    kept
geneX      geneX:3-26:132-154  CTGCTGTCGGACTCCTAGTTACG  GTCTTACGGACGGGAAAGAACC 151           1                  True    kept
...
```

`result/`에는 `candidates.tsv` 외에도 각 프라이머 쌍이 결합하는 위치를 담은
`amplicons.tsv`, 실행 정보를 담은 `manifest.json`, 병합된 설정 파일인 `run.yaml`이
생성됩니다. 실제 분석에서는 참조 서열을 본인의 유전체나 전사체로 바꾸고, 원하는 방식으로
타겟을 지정하면 됩니다.

## 사용법

아래 예제에서 `reference.fa`는 사용자가 준비한 FASTA입니다. 유전체, 전사체, 또는 다른
multi-FASTA를 사용할 수 있습니다. 모든 실행은 기본적으로 같은 구조를 가집니다.

```bash
ehtpcr run -r reference.fa  <타겟 지정 방법>  -o result/
```

옵션을 명령줄 플래그로 줄 수도 있고, YAML 설정 파일(`--config`)에 담을 수도 있습니다.
둘을 함께 쓰면 명령줄 플래그가 설정 파일 값을 덮어씁니다. 전체 옵션은 다음 명령으로 확인할
수 있습니다.

```bash
ehtpcr run --help
```

도움말은 reference, target, design, runtime 범주로 나뉘어 표시됩니다.

### 1. 타겟 지정

프라이머를 설계할 타겟은 네 가지 방식으로 지정할 수 있습니다. 가진 정보에 맞는 방식을
고르면 됩니다. 보통은 사용한 플래그에 따라 모드가 자동으로 결정되므로 `--target-mode`를
직접 지정할 일은 많지 않습니다. 기본 모드는 `--name`입니다.

선택 조건에 맞는 타겟이 **하나도 없으면** 실행은 에러로 멈춥니다. 보통은 오타나 잘못된
ID를 빨리 발견하는 편이 안전하기 때문입니다. 빈 결과가 정상일 수 있는 배치 파이프라인에서는
`--no-fail-on-no-targets`를 사용하세요.

**이름으로 지정**: 유전자나 전사체 ID를 알고 있을 때 사용합니다. 기본 모드입니다.
FASTA header에 지정한 문자열이 *포함된* 항목을 모두 가져오므로, 같은 유전자의 isoform도
함께 타겟이 됩니다.

```bash
ehtpcr run -r reference.fa --name Glyma.08G110400 -o result/
```

예를 들어 FASTA에 `>Glyma.08G110400.1 ...`과 `>Glyma.08G110400.2 ...`가 있으면 두 항목
모두 타겟으로 사용됩니다. 여러 유전자를 지정하려면 `--name`을 반복합니다.

```bash
ehtpcr run -r reference.fa --name Glyma.08G110400 --name Glyma.05G001200 -o result/
```

**좌표로 지정**: 타겟 위치를 알고 있을 때 사용합니다. `--region` 좌표는 genome browser처럼
1-based inclusive 형식입니다. `:+` 또는 `:-` strand 표기는 선택 사항입니다.

```bash
ehtpcr run -r genome.fa --region "Chr08:110340-112000:+" -o result/
```

**유전자 ID와 GFF로 지정**: GFF3 gene model이 있을 때 사용합니다. `--locus`와 `--gff`를
함께 주면 eHT-PCR이 해당 `ID=` feature를 찾아 그 구간을 타겟으로 사용합니다.

```bash
ehtpcr run -r genome.fa --gff genes.gff3 --locus AT1G01010 -o result/
```

**서열로 지정**: 타겟 이름이나 좌표는 모르고 서열만 있을 때 사용합니다. `--query-fasta`를
주면 eHT-PCR이 정렬을 통해 참조 서열 안의 위치를 찾습니다. 같은 쿼리가 두 곳 이상에
매치되면 실행을 중단합니다. 타겟 위치가 모호한 상태에서 프라이머를 설계하면 결과 해석이
위험하기 때문입니다.

```bash
ehtpcr run -r genome.fa --query-fasta my_gene.fa -o result/
```

쿼리는 DNA 또는 단백질 서열일 수 있습니다. 사용할 정렬기는 자동으로 선택됩니다. 자세한
내용은 [쿼리 서열 위치 찾기](#3-쿼리-서열-위치-찾기-target_fasta)를 참고하세요.

### 2. 특이성 검사와 조정

기본적으로 특이성은 설계에 사용한 참조 서열을 기준으로 검사합니다. 하지만 설계용 참조와
특이성 검사용 참조를 분리하고 싶은 경우가 많습니다. 예를 들어 기존 soybean workflow에서는
전사체에서 프라이머를 설계하고, 알려진 중복 유전자를 제거한 curated genome을 기준으로
특이성을 검사했습니다. 이럴 때는 `--spec-reference`를 사용합니다.

```bash
ehtpcr run -r Gmax.transcript.fa --name Glyma.08G110400 \
  --spec-reference Gmax_curated.fa --max-mismatch 1 -o result/
```

`--max-mismatch`는 특이성 검사에서 가장 중요한 조정값입니다. off-target 결합 위치를 찾을 때
`bwa`가 허용할 mismatch 개수를 뜻합니다. 원래 workflow에서는 촘촘한 유전자 cluster에는
`1`, 그 외에는 `5`를 사용했습니다. 타겟별로 값을 다르게 주고 싶다면 설정 파일을 쓰면
됩니다.

```yaml
# myrun.yaml
reference:
  fasta: Gmax.transcript.fa
target:
  mode: name
  name: [Glyma.08G110400, Glyma.05G001200]
specificity:
  reference:
    fasta: Gmax_curated.fa     # 별도의 curated reference, 선택 사항
  keep: unique                 # unique | all | nomatch
  params:
    bwa_aln:
      max_mismatch: 5
      overrides:               # 까다로운 cluster에는 더 엄격한 mismatch 기준 적용
        Glyma.08G110400.1: 1
runtime:
  threads: 8
  jobs: 4                      # 여러 타겟을 병렬 처리
```

```bash
ehtpcr run --config myrun.yaml -o result/
```

`keep`는 어떤 쌍을 `kept`로 출력에 남길지 결정합니다. 기본값 `unique`는 on-target 증폭
산물이 하나인 쌍만, `all`은 모든 후보를 판정과 함께, `nomatch`는 증폭 산물이 전혀 없는
쌍만 남깁니다.

특이성 검사를 건너뛰고 프라이머 후보만 만들고 싶다면 `--no-spec`를 추가하세요.

### 3. 쿼리 서열 위치 찾기 (`target_fasta`)

`--query-fasta` 또는 `target.mode: target_fasta`로 쿼리 서열을 주면, `--locator`가 그
서열을 참조 서열에 어떻게 정렬할지 결정합니다. 기본값 `auto`는 쿼리 종류를 판단해 DNA에는
`blastn`, 단백질에는 `tblastn`을 사용합니다.

한 파일 안에서는 하나의 locator만 선택됩니다. DNA와 단백질 레코드를 한 파일에 섞지 않는
것이 좋습니다. 꼭 섞어야 한다면 `--locator`를 명시하세요. 중복된 쿼리 이름은 허용되지
않습니다.

| `--locator` | 적합한 쿼리 | 설명 |
|-------------|-------------|------|
| `auto` | 모든 쿼리, 기본값 | DNA → `blastn`, 단백질 → `tblastn` |
| `bwa` | 참조 서열과 거의 동일한 DNA | 추가 의존성 없음 |
| `blastn` | 차이가 있을 수 있는 DNA | 민감도가 높고 paralog 구분에 적합함. 전사체-유전체 비교도 처리 가능 |
| `minimap2` | 유전체에 매핑할 spliced transcript | intron-aware. paralog 판별용으로는 **부적합** |
| `tblastn` | 유전체에 매핑할 단백질 | 6-frame translation 사용. paralog 구분에 적합함 |
| `miniprot` | 유전체에 매핑할 단백질 | splice/frameshift-aware. paralog 판별용으로는 **부적합** |

```bash
ehtpcr run -r genome.fa --query-fasta my_protein.fa --locator tblastn -o result/
```

paralog를 구분하는 것이 중요하다면 `blastn` 또는 `tblastn`을 사용하세요. 정확한 spliced
gene structure가 필요한 경우에만 `minimap2`나 `miniprot`을 선택하는 편이 좋습니다.

### 4. 설계 조건 조정

모든 후보 프라이머는 아래 조건을 통과해야 합니다. 기본값은 대부분의 genotyping / qPCR
설계에 적합하며, CLI 옵션이나 설정 파일의 `design:` 블록에서 바꿀 수 있습니다.

| 조건 | 플래그 | 기본값 | 의미 |
|------|--------|--------|------|
| 길이 | `--len-min` / `--len-max` | 18–24 nt | 프라이머 길이 범위 |
| GC 함량 | `--gc-min` / `--gc-max` | 0.30–0.70 | G/C 비율(너무 낮으면 결합이 약하고, 너무 높으면 변성이 어려움) |
| Tm | `--tm-min` / `--tm-max` | 57–62 °C | 녹는점 — [녹는점](#녹는점-melting-temperature) 참고 |
| 산물 크기 | `--product-min` / `--product-max` | 150–250 bp | 증폭 산물 길이 |
| 3' GC-clamp | `--gc-clamp` | 2 | 결합을 안정화하기 위해 3' 끝에 요구하는 G/C 염기 수(0이면 사용 안 함) |

또한 프라이머 쌍의 정방향·역방향 프라이머는 서로 **잘 맞아야** 합니다. 기본적으로 두
프라이머의 Tm 차이는 최대 1 °C(`design.pair.tm_diff.max`)까지 허용되며, GC 차이도 같은
방식으로 `design.pair.gc_diff.max`로 제한합니다. 균일한 프라이머 세트가 필요하면 이 값을
좁히고, 더 많은 후보 쌍이 필요하면 넓히세요.

## 출력

출력 디렉터리에는 다음 파일이 생성됩니다.

- **`candidates.tsv`**: 설계된 프라이머 쌍이 한 행씩 기록됩니다. 두 프라이머 서열, 좌표,
  Tm/GC, 생성된 유효 증폭 산물 수(`n_valid_amplicons`), `unique` 여부, `keep_decision`이
  포함됩니다. **실제로 사용할 후보는 `keep_decision == kept`인 행**입니다. 이 행은
  `unique == True`인 프라이머 쌍입니다.
- **`amplicons.tsv`**: 발견된 증폭 산물이 한 행씩 기록됩니다. 위치와 프라이머별 mismatch
  개수가 포함되며, `unique` 판정의 근거가 됩니다.
- **`manifest.json`**: eHT-PCR과 `bwa` 버전, 참조 서열 checksum, 실제 명령줄, timestamp,
  타겟별 mismatch override가 기록됩니다.
- **`run.yaml`**: 설정 파일과 명령줄 옵션을 병합한 최종 설정입니다. 실행을 재현하는 데
  필요한 정보를 담고 있습니다.

## 성능

가장 시간이 많이 드는 단계는 특이성 검사, 즉 `bwa` 단계입니다. 성능에 영향을 주는 주요
옵션은 두 가지입니다.

- `--threads`: 각 `bwa` 호출에 사용할 thread 수
- `--jobs`: 동시에 처리할 타겟 수

타겟이 많을 때는 `--jobs`가 특히 중요합니다. `--jobs` 값이 달라도 출력 결과는 동일합니다.

`bwa`, BLAST+, miniprot 인덱스와 FASTA `.fai`, GFF database는 한 번 만든 뒤 캐시됩니다.
캐시 키는 각 참조 파일의 경로, 크기, mtime을 기준으로 만들어지므로 같은 참조 서열을
반복해서 사용할 때는 기존 인덱스를 재사용합니다.

기본 캐시 위치는 `$XDG_CACHE_HOME/ehtpcr` 또는 `~/.cache/ehtpcr`입니다. 컨테이너나 공유
HPC 노드처럼 이 위치가 읽기 전용인 환경에서는 전체 캐시를 쓸 수 있는 위치로 바꾸세요.

```bash
ehtpcr run ... --cache-dir /scratch/you/ehtpcr
```

또는 환경 변수를 사용할 수 있습니다.

```bash
export EHTPCR_CACHE_DIR=/scratch/you/ehtpcr
```

캐시 위치에 쓸 수 없으면 실행은 `"cache directory is not writable"` 에러로 멈춥니다.

## 녹는점 (Melting temperature)

nearest-neighbor Tm은 염 농도와 프라이머 농도에 영향을 받습니다. eHT-PCR은 모델 간 값을
비교할 수 있도록 모든 모델에서 농도를 **50 mM Na+ / 200 nM**으로 고정합니다.

기본값인 `legacy` 모델은 기존 도구의 계산기(Breslauer/Borer parameter)를 그대로 옮긴
것입니다. 기존 parameter set이 같은 방식으로 동작하도록 기본값으로 유지했습니다. 새로
설계하는 경우에는 사실상 표준인 SantaLucia 1998 계산을 사용하는 **`primer3`** 또는
`santalucia2004`를 권장합니다.

그 밖에 `breslauer`, `sugimoto`, `santalucia`, `wallace`, `gc` 모델도 `--tm-model` 또는
`design.tm.model`로 선택할 수 있습니다. `legacy` 모델은 SantaLucia 계열 모델보다 Tm을 몇
°C 높게 계산하는 경향이 있으므로, 모델을 바꾸면 Tm 범위를 다시 확인하세요.

## 참고

**특이성 판정 방식.** eHT-PCR의 unique 판정은 열역학 계산이 아니라 *정렬 결과*를 기준으로
합니다. 각 프라이머는 `bwa aln`으로 매핑되며, `max_mismatch` 이하의 mismatch를 위치와
상관없이 허용합니다.

실제 PCR에서는 **3' 말단** mismatch가 특히 중요합니다. 3' 끝이 맞지 않는 프라이머는 대개
연장되지 않지만, 단순 mismatch 개수만으로는 이런 차이를 충분히 반영할 수 없습니다. 따라서
`amplicons.tsv`는 "이 프라이머들이 결합할 가능성이 있는 위치" 목록으로 해석하고,
타겟별로 `max_mismatch`를 조정하며, 가능하면 curated specificity reference를 사용해
설계하는 것이 좋습니다. 특이성 엔진은 교체 가능한 구조라서, 이후 열역학 기반 백엔드를
추가할 수 있습니다.

**`name` 모드의 on-target 범위.** `name` 모드에서는 matched FASTA entry 위에 생긴
증폭 산물을 on-target으로 간주합니다. 하위 좌표까지 따로 검사하지는 않습니다. 전사체에서는
전사체 하나가 분석 단위이므로 이 방식이 자연스럽습니다.

반면 유전체 FASTA에서 염색체 전체를 이름으로 선택하면, 같은 염색체의 다른 위치에서 생긴
산물도 on-target으로 처리될 수 있습니다. 유전체에서 contig 내부의 특정 구간을 타겟으로
삼고 싶다면 `region` 또는 `locus` 모드를 사용하세요.

**좌표.** 내부 좌표는 모두 0-based half-open 형식입니다. TSV 열 이름도 `start0`, `end0`처럼
이를 명시합니다. 반대로 `--region`과 GFF feature처럼 관례적으로 1-based를 쓰는 입력은
경계에서 변환됩니다. 사용자는 genome browser에서 보는 방식 그대로 좌표를 입력하면 됩니다.

## 라이선스

eHT-PCR 자체 코드는 MIT 라이선스입니다. 자세한 내용은 [LICENSE](LICENSE)를 참고하세요.
다만 핵심 의존성 중 하나인 [`primer3-py`](https://github.com/libnano/primer3-py)는
**GPLv2**입니다. `primer3` Tm 모델에서 사용되며, eHT-PCR을 설치하면 함께 설치됩니다.

완전히 허용적인(permissive) 의존성 구성이 필요하다면 `primer3` 모델을 사용하지 않으면 됩니다.
기본 `legacy` 모델, Biopython nearest-neighbor 모델, Wallace, GC 모델은 permissive
라이선스 패키지에만 의존합니다.
