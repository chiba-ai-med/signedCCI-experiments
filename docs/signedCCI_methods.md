---
title: "符号付き細胞間相互作用（signed cell–cell interaction）解析：材料と方法"
author: "signedCCI-experiments"
date: "2026-06-12"
documentclass: article
geometry: margin=2.4cm
mainfont: "Noto Sans CJK JP"
CJKmainfont: "Noto Sans CJK JP"
fontsize: 10pt
linestretch: 1.15
---

## 0. 本書の位置づけと用語に関する注意

本書は、リポジトリ `signedCCI-experiments` で行った一連の解析の「材料と方法
（Materials and Methods）」をまとめたものである。読者が用語でつまずかないよう、
**専門用語は初出時に定義し、巻末に用語解説（第11節）を置く**。

特に、本書で用いる用語は次の3種類に区別して扱う。誤解を避けるため、本文中でも
都度その区別を明示する。

- **確立した用語**：分野で一般的に使われる語（例：細胞間相互作用、リガンド–レセプター）。
- **ソフトウェアが定義する名称**：`symTensor` / `scTensor` パッケージの関数名が
  そのまま概念名になっているもの（例：`SignedPageRank`、`BuildSignedCCI`）。
  これらは「その関数が実装する具体的アルゴリズム」を指し、文献上の標準用語とは限らない。
- **本研究の記述的造語**：二つの概念を区別するために本研究で導入した言い回し
  （例：「機能アウトカム符号」と「分子シグナル符号」、リソース名「signedLRBase」）。
  造語であることを明記する。

---

## 1. 概要

単一細胞 RNA-seq データにおける細胞間相互作用（cell–cell interaction; 以下 CCI）を、
**符号付き（signed）**、すなわち各リガンド–レセプター（ligand–receptor; 以下 LR）ペアに
「活性化など正の作用（$+1$）」「抑制・細胞死誘導など負の作用（$-1$）」の符号を割り当てた
有向グラフとして表現し、その上で符号付きランダムウォーク（`SignedPageRank`）と
符号付き経路寄与（`SignedPathContribution`）を計算した。中心的な性質は
$(-)\times(-)=(+)$（「敵の敵は味方」）であり、二重の抑制を経由した**間接的な正の
相互作用**を検出できる。

検証は (i) 人工（玩具）データでの数理的性質の確認、(ii) 免疫×腫瘍の実データ
（メラノーマ腫瘍微小環境）への適用、(iii) 符号を機能カテゴリから割り当てる
**signedLRBase**（本研究の造語）の構築、(iv) 符号の文脈依存性の評価、の順で行った。

---

## 2. 計算環境とソフトウェア

すべての解析は `conda` 仮想環境 `signedcci`（R 4.5.3）で実行した。主要パッケージと
その入手元を表1に示す。`symTensor` と `scTensor` は Bioconductor リリース版ではなく、
開発元 GitHub（`rikenbit`）の開発版を用いた。両者の符号付き CCI 関数
（`BuildSignedCCI`、`SignedPageRank`、`SignedPathContribution`）は素の base R 実装で
あり、重い Bioconductor 依存は計算自体には不要である（インストール時の依存解決にのみ必要）。

表1：主要ソフトウェア

| パッケージ／資源 | バージョン | 入手元・役割 |
|---|---|---|
| R | 4.5.3 | conda-forge |
| symTensor | 0.1.0 | `rikenbit/symTensor`（GitHub, branch `main`）。`SignedPageRank` 等 |
| scTensor | 2.23.1 | `rikenbit/scTensor`（GitHub, branch `master`）。`BuildSignedCCI` |
| rTensor | 1.5.0 | テンソル演算（symTensor の依存） |
| progeny | 1.32.0 | Bioconductor。経路活性推定（第10節） |
| OmniPath | REST API | `omnipathdb.org`。符号照合（第7節） |

インストール手順（要約）：`mamba create -n signedcci -c conda-forge -c bioconda
bioconductor-sctensor=2.20.0 r-remotes r-biocmanager git` で重い依存ツリーを一括
解決した後、`remotes::install_github("rikenbit/symTensor", ref="main")`、続いて
`remotes::install_github("rikenbit/scTensor", ref="master")` で開発版を上書き導入した
（`scTensor` 2.23.1 は `symTensor` に依存するため、導入順序は symTensor → scTensor）。

---

## 3. データセット

実データには **GSE72056**（Tirosh ら 2016、転移性メラノーマの腫瘍微小環境
〔tumor micro-environment; TME〕の単一細胞 RNA-seq）を用いた。GEO の補足ファイル
`GSE72056_melanoma_single_cell_revised_v2.txt.gz` を取得した。発現値は
$\log_2(\mathrm{TPM}/10 + 1)$ 形式である（TPM = transcripts per million、転写産物存在量の
正規化単位）。

ファイル先頭の4行はアノテーション行であり、第1行が細胞 ID、第2行が腫瘍 ID、
第3行が悪性度コード（malignant：1=非悪性, 2=悪性, 0=判定不能）、第4行が
非悪性細胞型コード（1=T, 2=B, 3=Macrophage, 4=Endothelial, 5=CAF, 6=NK）である。
重複する遺伝子シンボルは行ごとの最大値で代表させて統合し、最終的に **23,684 遺伝子
× 4,645 細胞**の行列を得た。

---

## 4. 細胞型の割り当て

各細胞を次の規則で細胞型に割り当てた。

1. 悪性度コード = 2 の細胞を **Tumor** とする。
2. 悪性度コード = 1（非悪性）の細胞は、細胞型コードに従い
   **Bcell / Macrophage / Endo / CAF / NK** に展開する。
3. T 細胞（コード 1）はマーカー遺伝子の発現で細分する：`FOXP3` が発現（$>0$）なら
   **Treg**（制御性 T 細胞）、そうでなく `CD8A` または `CD8B` が発現なら **CD8**
   （細胞傷害性 T 細胞）、いずれでもなければ **CD4conv**（従来型 CD4 T 細胞）とする。

細胞数が 20 未満の細胞型は除外した。得られた 9 細胞型と細胞数を表2に示す。

表2：細胞型と細胞数

| 細胞型 | n | 細胞型 | n |
|---|---|---|---|
| Tumor | 1257 | Macrophage | 119 |
| CD8 | 1188 | Endo | 62 |
| CD4conv | 692 | CAF | 56 |
| Bcell | 512 | NK | 51 |
| Treg | 160 | | |

各細胞型について、構成細胞にわたる平均 $\log$ 発現を計算し、遺伝子 × 細胞型の
平均発現行列を構築した（この行列を以後 $E$ と書く）。この手続きはコード上では
共通関数 `R/gse72056.R` に切り出してある。

---

## 5. 符号付き通信強度の定義（BuildSignedCCI）

LR ペア $k$ について、送り手細胞型 $i$ のリガンド発現を $L_{ik}$、受け手細胞型 $j$ の
レセプター発現を $R_{jk}$ とし、ペアごとの通信強度を

$$ X_{ijk} = L_{ik}\, R_{jk} $$

で定義する（$i$ = sender、$j$ = receiver、$k$ = LR ペア）。各 LR ペアには符号
$\sigma_k \in \{-1, +1\}$ が付与される（第6–7節）。これを符号別に集約して、正の
通信行列 $A^{+}$ と負の通信行列 $A^{-}$ を得る：

$$ A^{+}_{ij} = \!\!\sum_{k:\,\sigma_k = +1}\!\! X_{ijk}, \qquad
   A^{-}_{ij} = \!\!\sum_{k:\,\sigma_k = -1}\!\! X_{ijk}. $$

すなわち LR 次元は集約され、以後の解析は $A^{+}, A^{-}$ のみを用いる。これは
`scTensor::BuildSignedCCI(lig_expr, rec_expr, lr_sign)` が実装する処理である
（出力：`A_pos`, `A_neg`, ペアごとの $X$ テンソル、各エッジの LR 内訳 `edge_table`）。

入力行列 `lig_expr` / `rec_expr` は（細胞型 × LR ペア）形式で、ヘルパー関数
`build_lr_matrices()` が構築する。人工データ（玩具）モードでは各 LR ペアの送り手・
受け手のみを発現させて既知の有向エッジを再現し、実データモードでは平均発現行列
$E$ から各 LR ペアのリガンド／レセプター遺伝子の発現を引く。

**向きに関する重要な注意（実装上の規約）：** `BuildSignedCCI` の出力 $A_{ij}$ は
$i$=送り手（行）・$j$=受け手（列）、すなわち $A[\text{from}, \text{to}]$ である。
一方、後述の `SignedPageRank` / `SignedPathContribution` は列方向の正規化
（$\mathrm{colSums}$）を用い、$A[\text{to}, \text{from}]$ を前提とする。したがって
両者を接続する際は **転置 $A^{\top}$ を渡す**。この処理は `signed_cci()` 内で行っている。

---

## 6. 符号付き PageRank（SignedPageRank）

**用語：** ここでいう「Signed PageRank」は `symTensor` パッケージの関数
`SignedPageRank` が実装する具体的アルゴリズムの名称であり、状態空間を正負に二重化した
符号付きランダムウォークである。標準的な PageRank（Web ページの重要度を、リンクを
たどるランダムウォークの定常分布で測る手法）を、符号付きグラフに拡張したものである。

出次数（out-degree、各ノードから出る重みの総和）を

$$ D_i = \sum_j \bigl( A^{+}_{ij} + A^{-}_{ij} \bigr) $$

とし、列正規化した遷移確率行列を $P^{+}_{ij} = A^{+}_{ij}/D_i$、
$P^{-}_{ij} = A^{-}_{ij}/D_i$ とする。正・負それぞれの状態をもつ拡張状態
$R = (r^{+}, r^{-})^{\top}$ を導入し、遷移行列を

$$ M = \begin{pmatrix} P^{+} & P^{-} \\ P^{-} & P^{+} \end{pmatrix} $$

と構成する。$M$ の非対角ブロックが $P^{-}$ である点が鍵で、負のエッジをたどると
正・負の状態が入れ替わる。これにより**負を2回たどると正に戻る**
（$(-)\times(-)=(+)$）。反復は

$$ R^{(t+1)} = \alpha\, M^{\top} R^{(t)} + (1-\alpha)\, V $$

で収束まで行う。$\alpha$ は減衰係数（damping factor、ランダムウォークがリンクを
たどる確率。残りの $1-\alpha$ で初期分布へテレポートする）で、本研究では
$\alpha = 0.85$ を用いた。$V$ は個人化ベクトル（personalization vector、テレポート
先の分布。既定では一様）。収束判定は L1 ノルムの変化が許容誤差 $\mathrm{tol}=10^{-8}$ を
下回ること、最大反復数 1000 とした。出次数 0 のノード（dangling node）の質量は
個人化ベクトルへ再配分する。

出力は各細胞型について 4 つの量：`positive`（正の定常スコア）、`negative`（負の
定常スコア）、`net` $=$ positive $-$ negative、`total` $=$ positive $+$ negative。

---

## 7. 符号付き経路寄与（SignedPathContribution）と寄与 LR の分解

**用語：** 「Signed Path Contribution」も `symTensor` の関数 `SignedPathContribution`
が実装する処理の名称であり、符号付きグラフ上の経路を列挙し、各経路について符号と
寄与量を返す。

ノード間の経路（hop 数 1〜`max_hop`、本研究では `max_hop=3`）を列挙し、各経路に
ついて以下を返す：始点 `source`、終点 `target`、hop 数、経路上のノード列 `path`、
各 hop の符号列 `signs`、経路全体の符号 `net_sign`（各 hop 符号の積）、寄与量
`contribution`（経路上の遷移確率の積）。遷移確率は第6節と同じ列正規化 $P^{\pm}$ を
用いる。中心的性質は、**負のエッジを偶数回含む経路は $net\_sign = +1$ となる**こと
である。例：$\mathrm{Treg}\xrightarrow{-}\mathrm{CD8}\xrightarrow{-}\mathrm{Tumor}$ は
$net\_sign = +1$。

**寄与 LR の分解：** `BuildSignedCCI` の `edge_table` は各エッジを構成する LR ペア
ごとの $X_{ijk}$ を保持する。これを用いて、各符号付きエッジ $(i\to j)$ を寄与の
大きい LR ペア順にランキングし、主要 LR を同定した
（`experiments/real/lr_contribution.R`）。

---

## 8. 符号付き LR データベース signedLRBase の構築

**用語（造語の明示）：** 「signedLRBase」は本研究で構築した符号付き LR リソースの
名称（造語）である。また本研究では符号の意味を二つに区別し、それぞれを記述的造語で
呼ぶ：

- **機能アウトカム符号（functional-outcome sign）**〔本研究の造語〕：受け手細胞に
  とっての作用の向き。$+1$ = 活性化・生存促進・増殖促進、$-1$ = 免疫抑制・疲弊誘導・
  細胞死誘導。**本研究が必要とするのはこちら**。
- **分子シグナル符号（molecular-signalling sign）**〔本研究の造語〕：リガンドが
  レセプターの**シグナル伝達**を活性化（agonist）するか阻害（antagonist）するか。
  既存の符号付き相互作用 DB（SIGNOR、OmniPath の `is_stimulation`/`is_inhibition`、
  CellPhoneDB v5）が持つのはこちら。

両者は一致しない。代表例：$\mathrm{FASLG}\to\mathrm{FAS}$ は分子的には agonist
（$+$）だが、細胞アウトカムは死（$-$）。したがって既存 DB の分子符号をそのまま
機能符号に流用できない。これが signedLRBase を構築する根拠である。

**構築手順（ルール＋DB照合）：**

1. **キュレーション（ルール層）**：免疫／TME に関連する高信頼 LR ペア 91 件を、機能
   カテゴリに基づいて手作業で符号付けした（`signed_lr_curated.csv`）。負（$-1$）の
   カテゴリ：免疫チェックポイント阻害性、NK 抑制性、細胞死、抑制性サイトカイン。
   正（$+1$）のカテゴリ：共刺激、活性化サイトカイン、ケモカイン（エフェクター
   動員）、増殖・生存因子、NK 活性化。文脈依存のペアには `context_dependent` フラグを
   付した。
2. **DB 照合層**：OmniPath の REST API（データセット `omnipath`,`ligrecextra`）から
   分子シグナル符号（`is_stimulation` / `is_inhibition`）を取得し、各ペアの
   ルール符号と照合した。
3. **出力**：`signedLRBase.csv`（列：ligand, receptor, sign, category,
   context_dependent, db_sign, agreement, confidence, source, notes）。最終符号は
   機能アウトカム符号（ルール）を採用し、分子符号との一致状況を併記する。

**照合結果：** 91 ペア中、分子符号と **一致 50**、**想定内 override 26**、
**DB 曖昧 15**、**真の要再検 0**。想定内 override とは、抑制性経路が「抑制性
レセプターを分子的に活性化する」ことで成り立つために生じる、機能 $-1$ vs 分子
$+1$ のコンフリクト（PD-L1/PD-1、CTLA4、TIGIT、HLA/NKG2A、TRAIL、TGFβ、IL10 等）で
ある。すなわち分子符号をそのまま使うと負ペアの約 29% を誤符号化する。既知アンカー
（PD-L1/PD-1 $=-1$、CD80/CD28 $=+1$、CD80/CTLA4 $=-1$、FASLG/FAS $=-1$ 等）で
符号を検証した。

実データへの適用では、この signedLRBase（91 ペア中 88 が当該データで使用可能）を
`build_lr_matrices()` に与えて第5–7節のパイプラインを実行した
（`experiments/real/run_real_signedLRBase.R`）。

---

## 9. 符号の文脈依存性の評価

機能アウトカム符号は受け手細胞型・状態・組織・濃度・時間に依存する関数であり、単一の
グローバル定数ではない。したがって signedLRBase は**事前分布（prior）**として扱い、
全文脈の網羅は目的としない。文脈依存性を二通りで評価した。

**(a) マーカーに基づく受け手状態チェック（ルールベース、一般化しない）。**
context_dependent ペアの受け手細胞型について、活性化マーカー群と疲弊／抑制マーカー群
（あるいは Macrophage の M1／M2 マーカー群）の発現を細胞型間で z 化して平均し、
活性化スコアと抑制スコアを比較した（`experiments/real/sign_refinement_markers.R`）。
**注意：** これは受け手の*全体状態*を測るものであり、特定リガンドの*因果作用*では
ない。またマーカー群は手作業のため、新しいデータセットごとに作り直しが必要で一般化
しない。本データでは CD8 が疲弊状態（疲弊スコア 1.57 > 活性化スコア 1.23）であった
（Tirosh らの古典的所見と一致）。

**(b) 作り置きモデル PROGENy（任意データに一般化）。**
PROGENy（**用語**：Perturbation-Response Genes の略。多数の摂動実験から学習された、
経路ごとの符号付き「応答遺伝子重み」を用いて経路活性を推定する確立した手法。出力は
経路の「フットプリント〔下流応答の足跡〕」活性）を、平均発現行列 $E$ に適用した
（`organism="Human"`, `top=100`, `scale=TRUE`、14 経路）。モデルは事前構築済みで
データセット固有の調整を要さないため、任意のヒト発現データにそのまま適用できる。
本データでは、TGFb 経路が CAF、EGFR 経路が Tumor、VEGF・JAK-STAT・NFkB 経路が
Macrophage、Trail 経路が B/NK で高く、生物学的に妥当に局在した。

**手続きの三段階と、消えない事前知識：** 任意データで符号を読む手続きは
(1) どの遺伝子を読むか、(2) 活性を測る、(3) 活性を $\pm$ 符号へ翻訳する、の3段から
なる。作り置きモデルを用いれば (1)(2) はデータ非依存だが、(3)（例：「TGFb 経路が
活性 $\Rightarrow$ 抑制 $-$」）には必ず薄い事前知識が残る。「活性 $\ne$ 細胞に
とって良い」ためで、これは FASLG が分子 agonist でも機能は死であることと同じ構図で
ある。なお NicheNet は本目的に使えない：そのリガンド–標的制御ポテンシャルは
**無符号**であり、活性化と抑制が加算的に寄与して向きを分離できないためである。

---

## 10. 検証（玩具データによるゲートテスト）

数理的性質は人工ネットワークで先に確認した（`experiments/synthetic/run_synthetic.R`）。
5 細胞型（Tumor, Treg, CD8, NK, DC）と生物学ラベル付き符号付き LR 表から、既知の
有向符号付きエッジを再現する `lig_expr`/`rec_expr` を構成し、以下を `stopifnot` で
検証した（すべて合格）：

- 2-hop 経路 $\mathrm{Treg}\xrightarrow{-}\mathrm{CD8}\xrightarrow{-}\mathrm{Tumor}$ の
  $net\_sign = +1$（$(-)\times(-)=(+)$）。
- 3-hop 経路
  $\mathrm{Tumor}\xrightarrow{+}\mathrm{Treg}\xrightarrow{-}\mathrm{CD8}\xrightarrow{-}\mathrm{Tumor}$
  の $net\_sign = +1$（自己強化的な抗腫瘍ループ）。
- 単一 $(-)$ ホップは $net\_sign = -1$、単一 $(+)$ ホップは $+1$。
- `A_pos`/`A_neg` が設計した有向符号付きエッジと一致、`SignedPageRank` の収束。

---

## 11. 用語解説

- **細胞間相互作用（CCI）**〔確立〕：細胞が分泌・膜タンパク質を介して互いに信号を
  やり取りすること。
- **リガンド–レセプター（LR）ペア**〔確立〕：信号を送る分子（リガンド）と受ける分子
  （レセプター）の組。
- **符号付き CCI**〔本研究の枠組み〕：各 LR ペアに正負の符号を付けて表現した CCI。
- **機能アウトカム符号 / 分子シグナル符号**〔本研究の造語〕：第8節参照。受け手細胞への
  作用の向き／受容体シグナルへの活性化・阻害の向き。両者は一致しないことがある。
- **`BuildSignedCCI` / `SignedPageRank` / `SignedPathContribution`**〔ソフトウェア
  定義〕：`scTensor`/`symTensor` の関数名。それぞれ符号付き通信行列の構築、符号付き
  ランダムウォーク、符号付き経路寄与を計算する。
- **signedLRBase**〔本研究の造語〕：機能アウトカム符号を付与した LR リソース（第8節）。
- **PageRank / 減衰係数 / 個人化ベクトル / dangling node**〔確立〕：第6節参照。
- **out-degree（出次数）**〔確立〕：あるノードから出るエッジ重みの総和。
- **net_sign / contribution**〔ソフトウェア定義〕：経路の符号（hop 符号の積）／
  経路上の遷移確率の積（第7節）。
- **TME（腫瘍微小環境）**〔確立〕：腫瘍細胞と免疫・間質細胞などが共在する局所環境。
- **Treg / CD8 / CAF / NK / Endo**〔確立〕：制御性 T 細胞／細胞傷害性 T 細胞／がん
  関連線維芽細胞／ナチュラルキラー細胞／内皮細胞。
- **フットプリント（footprint）**〔確立〕：上流のシグナル活性を、その下流で変動する
  遺伝子群の発現パターンから推定する考え方。
- **PROGENy / DoRothEA / CytoSig**〔確立〕：いずれも作り置きの参照（経路応答遺伝子／
  符号付き転写因子制御セット／サイトカイン応答シグネチャ）を用いて、任意データから
  活性を推定するツール。
- **事前分布・事前知識（prior）**〔確立〕：データを見る前に持っている知識。本書では
  signedLRBase の既定符号を「事前分布」として扱う、の意で用いる。

---

## 12. データ・コードの可用性

コードと結果は `signedCCI-experiments`（GitHub: `chiba-ai-med/signedCCI-experiments`）
に置いた。主要スクリプト：`R/signed_lr_utils.R`（パイプライン中核）、
`R/signed_lrbase.R`（signedLRBase ローダ）、`R/gse72056.R`（細胞型割り当て）、
`experiments/synthetic/run_synthetic.R`（玩具ゲートテスト）、
`experiments/real/run_real.R` および `run_real_signedLRBase.R`（実データ）、
`experiments/real/lr_contribution.R`（寄与 LR 分解）、
`experiments/signedLRBase/build_signedLRBase.R`（signedLRBase 構築）、
`experiments/real/sign_refinement_markers.R` / `progeny_demo.R`（文脈依存の評価）。
GSE72056 の生データは GEO から取得し、リポジトリには含めない（大容量のため）。

---

## 13. 参考文献

以下は標準的な出典である（正式投稿前には書誌情報の最終確認を推奨する）。

1. Tirosh I, *et al.* Dissecting the multicellular ecosystem of metastatic melanoma
   by single-cell RNA-seq. *Science*. 2016;352(6282):189–196. (GSE72056)
2. Tsuyuzaki K, *et al.* scTensor: Detection of cell–cell interaction from
   single-cell RNA-seq dataset by tensor decomposition. Bioconductor package；
   関連プレプリント bioRxiv 566182 (2019).
3. Türei D, *et al.* Integrated intra- and intercellular signaling knowledge for
   multicellular omics analysis. *Mol Syst Biol*. 2021;17(3):e9923. (OmniPath)
4. Türei D, Korcsmáros T, Saez-Rodriguez J. OmniPath: guidelines and gateway for
   literature-curated signaling pathway resources. *Nat Methods*. 2016;13:966–967.
5. Schubert M, *et al.* Perturbation-response genes reveal signaling footprints in
   cancer gene expression. *Nat Commun*. 2018;9:20. (PROGENy)
6. Browaeys R, Saeys Y, *et al.* NicheNet: modeling intercellular communication by
   linking ligands to target genes. *Nat Methods*. 2020;17:159–162.
7. Garcia-Alonso L, *et al.* Benchmark and integration of resources for the
   estimation of human transcription factor activities. *Genome Res*.
   2019;29:1363–1375. (DoRothEA)
8. Jiang P, *et al.* Systematic investigation of cytokine signaling activity at the
   tissue and single-cell levels. *Nat Methods*. 2021;18:1181–1191. (CytoSig)
9. Efremova M, *et al.* CellPhoneDB: inferring cell–cell communication from combined
   expression of multi-subunit ligand–receptor complexes. *Nat Protoc*.
   2020;15:1484–1506.
