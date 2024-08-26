document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_neighbor_result = await fetchAPI(
    "/api/preloadneighbor",
    0,
    csrftoken
  );
  console.log(preload_neighbor_result);
  document.getElementById("preloadresult").textContent =
    preload_neighbor_result.preload_result;
});
function select(s_ul, selected) {
  document.getElementById(s_ul + "_select").textContent = selected;
  document.getElementById(s_ul + "_input").value = selected;
}
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");

  const n_input = document.getElementById("n_input").value;
  const k_neighbor = document.getElementById("k_neighbor").value;
  const n_neighbor = document.getElementById("n_neighbor").value;
  data = {
    n_input: n_input,
    k_neighbor: k_neighbor,
    n_neighbor: n_neighbor,
  };

  const n_neighbor_rslt = await fetchAPI(
    "/api/neighbor",
    JSON.stringify(data),
    csrftoken
  );
  console.log(n_neighbor_rslt);
  toggleLoading(false, "processbutton");

  document.getElementById("imgbox").classList.remove("hidden");
  loadImage(
    "neighbor_result",
    "neighborhood_lmplot.png",
    "neighborhood_lmplot_container"
  );
  loadImage(
    "neighbor_result",
    "neighborhood_heatmap.png",
    "neighborhood_heatmap_container"
  );
}
