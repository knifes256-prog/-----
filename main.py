import json
from urllib.parse import quote

import requests
import streamlit as st
import streamlit.components.v1 as components

from redox_data import HALF_REACTIONS, SUBSTANCES, combine, parse, reactions_for, resolve

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
                "`redox_data.py`의 `HALF_REACTIONS`와 `SUBSTANCES`에 추가하면 확장할 수 있습니다."
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
