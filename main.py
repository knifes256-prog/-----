"""산화·환원 반쪽 반응 데이터 (25°C, 표준 환원 전위 E°, 단위 V).

형식: "산화형 + n e⁻ → 환원형"  (반응물 쪽 전자 계수가 이동 전자 수 n)
PubChem은 구조 정보만 제공하고 산화·환원 반응 데이터는 없으므로,
이 앱은 교과서에 나오는 대표 반쪽 반응을 직접 모아 둔 표를 사용합니다.
"""
import re
from collections import Counter
from math import gcd

HALF_REACTIONS = [
    ("F₂ + 2e⁻ → 2F⁻", 2.87),
    ("H₂O₂ + 2H⁺ + 2e⁻ → 2H₂O", 1.78),
    ("MnO₄⁻ + 8H⁺ + 5e⁻ → Mn²⁺ + 4H₂O", 1.51),
    ("Au³⁺ + 3e⁻ → Au", 1.50),
    ("Cl₂ + 2e⁻ → 2Cl⁻", 1.36),
    ("Cr₂O₇²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H₂O", 1.33),
    ("O₂ + 4H⁺ + 4e⁻ → 2H₂O", 1.23),
    ("Br₂ + 2e⁻ → 2Br⁻", 1.07),
    ("NO₃⁻ + 4H⁺ + 3e⁻ → NO + 2H₂O", 0.96),
    ("Ag⁺ + e⁻ → Ag", 0.80),
    ("Fe³⁺ + e⁻ → Fe²⁺", 0.77),
    ("I₂ + 2e⁻ → 2I⁻", 0.54),
    ("Cu²⁺ + 2e⁻ → Cu", 0.34),
    ("Sn⁴⁺ + 2e⁻ → Sn²⁺", 0.15),
    ("2H⁺ + 2e⁻ → H₂", 0.00),
    ("Pb²⁺ + 2e⁻ → Pb", -0.13),
    ("Ni²⁺ + 2e⁻ → Ni", -0.26),
    ("Fe²⁺ + 2e⁻ → Fe", -0.44),
    ("Zn²⁺ + 2e⁻ → Zn", -0.76),
    ("2H₂O + 2e⁻ → H₂ + 2OH⁻", -0.83),
    ("Al³⁺ + 3e⁻ → Al", -1.66),
    ("Mg²⁺ + 2e⁻ → Mg", -2.37),
    ("Na⁺ + e⁻ → Na", -2.71),
    ("Ca²⁺ + 2e⁻ → Ca", -2.87),
    ("K⁺ + e⁻ → K", -2.93),
    ("Li⁺ + e⁻ → Li", -3.04),
]

# 영어 이름(PubChem 검색용) → 한글 별칭, 이 앱이 아는 화학종
SUBSTANCES = {
    "iron": {"ko": ["철"], "species": ["Fe", "Fe²⁺", "Fe³⁺"]},
    "copper": {"ko": ["구리"], "species": ["Cu", "Cu²⁺"]},
    "zinc": {"ko": ["아연"], "species": ["Zn", "Zn²⁺"]},
    "silver": {"ko": ["은"], "species": ["Ag", "Ag⁺"]},
    "gold": {"ko": ["금"], "species": ["Au", "Au³⁺"]},
    "sodium": {"ko": ["나트륨", "소듐"], "species": ["Na", "Na⁺"]},
    "potassium": {"ko": ["칼륨", "포타슘"], "species": ["K", "K⁺"]},
    "lithium": {"ko": ["리튬"], "species": ["Li", "Li⁺"]},
    "calcium": {"ko": ["칼슘"], "species": ["Ca", "Ca²⁺"]},
    "magnesium": {"ko": ["마그네슘"], "species": ["Mg", "Mg²⁺"]},
    "aluminum": {"ko": ["알루미늄"], "species": ["Al", "Al³⁺"]},
    "nickel": {"ko": ["니켈"], "species": ["Ni", "Ni²⁺"]},
    "lead": {"ko": ["납"], "species": ["Pb", "Pb²⁺"]},
    "tin": {"ko": ["주석"], "species": ["Sn²⁺", "Sn⁴⁺"]},
    "hydrogen": {"ko": ["수소"], "species": ["H₂", "H⁺"]},
    "oxygen": {"ko": ["산소"], "species": ["O₂"]},
    "water": {"ko": ["물"], "species": ["H₂O"]},
    "hydrogen peroxide": {"ko": ["과산화수소"], "species": ["H₂O₂"]},
    "chlorine": {"ko": ["염소"], "species": ["Cl₂", "Cl⁻"]},
    "bromine": {"ko": ["브로민", "브롬"], "species": ["Br₂", "Br⁻"]},
    "iodine": {"ko": ["아이오딘", "요오드"], "species": ["I₂", "I⁻"]},
    "fluorine": {"ko": ["플루오린", "불소"], "species": ["F₂", "F⁻"]},
    "potassium permanganate": {"ko": ["과망가니즈산 칼륨", "과망간산칼륨"], "species": ["MnO₄⁻"]},
    "potassium dichromate": {"ko": ["다이크로뮴산 칼륨", "중크롬산칼륨"], "species": ["Cr₂O₇²⁻"]},
    "nitric acid": {"ko": ["질산"], "species": ["NO₃⁻"]},
}


def resolve(query: str):
    """입력(한글/영어)을 (PubChem 검색명, 내장 데이터 키 또는 None)으로 변환."""
    q = query.strip().lower()
    for en, info in SUBSTANCES.items():
        if q == en or q in [k.lower().replace(" ", "") for k in info["ko"]] or q in info["ko"]:
            return en, en
        if q.replace(" ", "") in [k.replace(" ", "") for k in info["ko"]]:
            return en, en
    return query.strip(), None


def _parse_side(text):
    items = []
    for term in text.split(" + "):
        term = term.strip()
        m = re.match(r"^(\d+)?(.+)$", term)
        items.append((int(m.group(1) or 1), m.group(2)))
    return items


def parse(rxn: str):
    """'a + n e⁻ → b' → (왼쪽 [(계수,종)], 오른쪽 [(계수,종)], n)"""
    left, right = rxn.split(" → ")
    left, right = _parse_side(left), _parse_side(right)
    n = next(c for c, s in left if s == "e⁻")
    left = [(c, s) for c, s in left if s != "e⁻"]
    return left, right, n


def fmt(side):
    return " + ".join((str(c) if c != 1 else "") + s for c, s in side)


def reactions_for(substance_key):
    """해당 물질의 화학종이 들어 있는 반쪽 반응 목록과 역할(산화제/환원제)."""
    species = SUBSTANCES[substance_key]["species"]
    out = []
    for rxn, e in HALF_REACTIONS:
        left, right, n = parse(rxn)
        roles = []
        for sp in species:
            if any(s == sp for _, s in left):
                roles.append((sp, "산화제 (환원됨)"))
            if any(s == sp for _, s in right):
                roles.append((sp, "환원제 (산화됨)"))
        if roles:
            out.append({"rxn": rxn, "E": e, "n": n, "roles": roles})
    return out


def combine(rxn_a: str, e_a: float, rxn_b: str, e_b: float):
    """두 반쪽 반응을 합쳐 전체 반응을 만든다. E°가 큰 쪽이 환원(양극), 작은 쪽이 산화(음극)."""
    if e_a < e_b:
        rxn_a, e_a, rxn_b, e_b = rxn_b, e_b, rxn_a, e_a
    la, ra, na = parse(rxn_a)          # 환원: la → ra
    lb, rb, nb = parse(rxn_b)          # 산화: rb → lb (역반응)
    L = na * nb // gcd(na, nb)
    ma, mb = L // na, L // nb
    left, right = Counter(), Counter()
    for c, s in la:
        left[s] += c * ma
    for c, s in ra:
        right[s] += c * ma
    for c, s in rb:
        left[s] += c * mb
    for c, s in lb:
        right[s] += c * mb
    for s in set(left) & set(right):   # H⁺, H₂O 등 양쪽에 있는 종 상쇄
        k = min(left[s], right[s])
        left[s] -= k
        right[s] -= k
    L_ = [(c, s) for s, c in left.items() if c > 0]
    R_ = [(c, s) for s, c in right.items() if c > 0]
    e_cell = e_a - e_b
    dG = -L * 96485 * e_cell / 1000    # kJ/mol
    return {
        "equation": f"{fmt(L_)} → {fmt(R_)}",
        "reduction": rxn_a, "oxidation": rxn_b,
        "n": L, "E_cell": e_cell, "dG": dG,
        "spontaneous": e_cell > 0,
    }


import json
from urllib.parse import quote

import requests
import streamlit as st
import streamlit.components.v1 as components


PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

st.set_page_config(page_title="물질 정보 · 산화환원", page_icon="🧪", layout="centered")


@st.cache_data(ttl=86400, show_spinner=False)
def get_properties(name: str):
    url = f"{PUG}/compound/name/{quote(name)}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
    r = requests.get(url, timeout=15)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()["PropertyTable"]["Properties"][0]


@st.cache_data(ttl=86400, show_spinner=False)
def get_sdf(cid: int):
    """3D 좌표 SDF를 우선 받고, 없으면(단일 원자 등) 2D SDF로 대체."""
    for rtype in ("3d", "2d"):
        r = requests.get(f"{PUG}/compound/cid/{cid}/record/SDF?record_type={rtype}", timeout=15)
        if r.status_code == 200:
            return r.text, rtype
    return None, None


def viewer_3d(sdf: str, height: int = 420):
    html = f"""
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <div id="v" style="width:100%;height:{height}px;position:relative;"></div>
    <script>
      const v = $3Dmol.createViewer("v", {{backgroundColor: "white"}});
      v.addModel({json.dumps(sdf)}, "sdf");
      v.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.28}}}});
      v.zoomTo(); v.spin(true); v.render();
    </script>"""
    components.html(html, height=height + 10)


st.title("🧪 물질 정보 · 분자구조 · 산화환원")
st.caption("물질 이름(한글/영어)을 입력하면 분자구조와 산화·환원 반응을 보여줍니다. 구조 데이터: PubChem (미국 NIH)")

query = st.text_input("물질 이름", placeholder="예: 철, 과산화수소, caffeine, potassium permanganate")

if query:
    search_name, key = resolve(query)

    with st.spinner("PubChem에서 불러오는 중..."):
        try:
            props = get_properties(search_name)
        except requests.RequestException as e:
            st.error(f"PubChem 요청 중 오류가 났습니다: {e}")
            st.stop()

    if props is None:
        st.warning("PubChem에서 찾지 못했습니다. 영어 이름으로 다시 입력해 보세요. (한글은 내장 목록의 물질만 지원)")
        st.stop()

    cid = props["CID"]
    st.subheader(props.get("IUPACName") or search_name)
    c1, c2, c3 = st.columns(3)
    c1.metric("분자식", props["MolecularFormula"])
    c2.metric("분자량", props["MolecularWeight"])
    c3.metric("PubChem CID", cid)

    tab_struct, tab_redox = st.tabs(["🔬 분자 구조", "⚡ 산화·환원 반응"])

    with tab_struct:
        left, right = st.columns(2)
        with left:
            st.markdown("**2D 구조**")
            st.image(f"{PUG}/compound/cid/{cid}/PNG", use_container_width=True)
        with right:
            sdf, rtype = get_sdf(cid)
            st.markdown("**3D 구조** (드래그로 회전)" if rtype == "3d" else "**구조** (3D 좌표 없음 → 2D)")
            if sdf:
                viewer_3d(sdf)
            else:
                st.info("구조 파일을 받지 못했습니다.")

    with tab_redox:
        if key is None:
            st.info(
                "이 물질은 내장 산화·환원 표에 없습니다. PubChem은 구조 정보만 제공하고 "
                "산화·환원 반응 데이터는 없기 때문에, 이 앱은 대표 물질의 반쪽 반응만 직접 모아 두었습니다. "
                "이 파일 위쪽의 `HALF_REACTIONS`와 `SUBSTANCES`에 추가하면 확장할 수 있습니다."
            )
            st.caption("지원 물질: " + ", ".join(f"{v['ko'][0]}({k})" for k, v in SUBSTANCES.items()))
        else:
            rows = reactions_for(key)
            st.markdown("#### 표준 환원 전위 (25°C, 1 M, 1 atm)")
            for row in rows:
                who = " / ".join(f"{sp}: {role}" for sp, role in row["roles"])
                st.markdown(f"- `{row['rxn']}`  \n  E° = **{row['E']:+.2f} V** · {who}")

            st.markdown("#### 다른 반쪽 반응과 합쳐 보기")
            options = [f"{r} (E° {e:+.2f} V)" for r, e in HALF_REACTIONS]
            mine = st.selectbox("이 물질의 반쪽 반응", [f"{r['rxn']} (E° {r['E']:+.2f} V)" for r in rows])
            partner = st.selectbox("짝이 될 반쪽 반응", options)
            ra, ea = next((r, e) for r, e in HALF_REACTIONS if f"{r} (E° {e:+.2f} V)" == mine)
            rb, eb = next((r, e) for r, e in HALF_REACTIONS if f"{r} (E° {e:+.2f} V)" == partner)
            if ra == rb:
                st.info("서로 다른 반쪽 반응을 골라 주세요.")
            else:
                res = combine(ra, ea, rb, eb)
                st.success(f"**전체 반응**: {res['equation']}")
                st.markdown(
                    f"- 환원(E° 큰 쪽): `{res['reduction']}`\n"
                    f"- 산화(E° 작은 쪽, 역반응): `{res['oxidation']}`\n"
                    f"- 이동 전자 수 n = **{res['n']}**\n"
                    f"- E°cell = **{res['E_cell']:+.2f} V**, ΔG° = **{res['dG']:.0f} kJ/mol**"
                )
                st.caption("E°cell > 0 이면 표준 상태에서 자발적." + (" → 자발적 반응입니다." if res["spontaneous"] else ""))

st.divider()
st.caption("구조: PubChem 제공 · 3D 렌더링: 3Dmol.js · 전위 값은 교과서 표준값(참고용)")
