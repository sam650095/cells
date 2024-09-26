let cluster_lst = [];
let cluster_lst_L = [];
document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_spatialanalysis_results = await fetchAPI(
    "/api/preloadspatialanalysis",
    0,
    csrftoken
  );
  document.getElementById("preloadresult").textContent =
    preload_spatialanalysis_results.data.preload_result;
  console.log(preload_spatialanalysis_results);
});
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const spatialanalysis_results = await fetchAPI(
    "/api/spatialanalysis",
    0,
    csrftoken
  );
  console.log(spatialanalysis_results);
  cluster_lst = spatialanalysis_results.data["cluster_list"];
  cluster_lst_L = spatialanalysis_results.data["cluster_list_L"];
  toggleLoading(false, "processbutton");
  document.getElementById("imgbox").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");

  // show dropdown menu
  adddp(
    spatialanalysis_results.data["columns_list"],
    spatialanalysis_results.data["cluster_list"],
    spatialanalysis_results.data["method_list"]
  );
  // show image
  addimg(spatialanalysis_results.data["filename"]);
}

async function addbtn(m) {
  const csrftoken = getCookie("csrftoken");

  let data;
  let container;
  if (m == "DH") {
    const dh_ul_input = document.getElementById("dh_ul_input").value;
    container = "distances-heatmap-container";
    data = {
      dh_ul_input: dh_ul_input,
    };
  } else if (m == "NP") {
    const np_ul_input = document.getElementById("np_ul_input").value;
    const npd_ul_input = document.getElementById("npd_ul_input").value;
    container = "numeric-plot-container";
    data = {
      np_ul_input: np_ul_input,
      npd_ul_input: npd_ul_input,
    };
  } else if (m == "IH") {
    const ih_ul_input = document.getElementById("ih_ul_input").value;
    const ihm_ul_input = document.getElementById("ihm_ul_input").value;
    container = "interactions-heatmap-container";
    data = {
      ih_ul_input: ih_ul_input,
      ihm_ul_input: ihm_ul_input,
    };
  } else if (m == "VP") {
    const vp_ul_input = document.getElementById("vp_ul_input").value;
    container = "voronoi-plot-container";
    data = {
      vp_ul_input: vp_ul_input,
    };
  }
  toggleLoading(true, m + "_btn");
  const rslt = await fetchAPI(
    `/api/addspatial/${m}`,
    JSON.stringify(data),
    csrftoken
  );
  console.log(rslt);
  toggleLoading(false, m + "_btn");
  loadImage("spatial_result", rslt["filename"], container, false);
}
