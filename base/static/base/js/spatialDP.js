function adddp(columns_list, cluster_list, method_list) {
  dp("dh_ul", columns_list);
  dp("np_ul", columns_list);
  dp("npd_ul", cluster_list);
  dp("ih_ul", columns_list);
  dp("ihm_ul", method_list);
  dp("vp_ul", columns_list);
}
function addimg(filename) {
  for (let i = 0; i < filename[0].length; i++) {
    loadImage("spatial_result", filename[0][i], "distances-heatmap-container");
  }
  loadImage("spatial_result", filename[1], "numeric-plot-container");
  loadImage("spatial_result", filename[2], "interactions-heatmap-container");
  loadImage("spatial_result", filename[3], "voronoi-plot-container");
}
function dp(s_ul, data_list) {
  const ul = document.getElementById(s_ul);
  ul.innerHTML = "";
  for (let i = 0; i < data_list.length; i++) {
    const li = document.createElement("li");
    li.classList.add("cursor-pointer");
    li.innerHTML = `
            <a onclick="select('${s_ul}', '${data_list[i]}')"
            class="block px-4 py-2 hover:bg-gray-100 m-2"
        >${data_list[i]}</a>`;
    ul.appendChild(li);
  }
}
function select(s_ul, selected) {
  document.getElementById(s_ul + "_select").textContent = selected;
  document.getElementById(s_ul + "_input").value = selected;
}
