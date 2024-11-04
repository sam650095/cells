document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_clustering_results = await fetchAPI(
    "/api/preloadclustering",
    0,
    csrftoken
  );
  show_clustering_method(preload_clustering_results.data);
  // grab step
  const grabstep_rslt = await grabsteps(`/getSteps/clustering/process/`);
  const grabstep_umap_add_image_rslt = await grabsteps(
    `/getSteps/clustering/umap_add_image/`
  );
  const grabstep_umap_add_l_image_rslt = await grabsteps(
    `/getSteps/clustering/umap_add_l_image/`
  );
  console.log(grabstep_umap_add_image_rslt);
  // clustering process
  if (grabstep_rslt.message != "notfound") {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    document.getElementById("nextbtn").classList.remove("hidden");
    // clustering process
    document.getElementById("preloadul_input").value =
      grabstep_rslt.input_values["chosen_method"];
    document.getElementById("preloadul_select").textContent =
      grabstep_rslt.input_values["chosen_method"];
    document.getElementById("n_neighbors").value =
      grabstep_rslt.input_values["n_neighbors"];
    document.getElementById("resolution").value =
      grabstep_rslt.input_values["resolution"];

    show_result();

    if (Object.keys(grabstep_umap_add_l_image_rslt.input_values).length != 0) {
      loadImage(
        "umap_cluster",
        "clustering_leidens_bysample.png",
        "addleidens-container"
      );
      document.getElementById("selectedbox").textContent =
        grabstep_umap_add_image_rslt.input_values["chosen_method"];
    }
    if (Object.keys(grabstep_umap_add_image_rslt.input_values).length != 0) {
      loadImage(
        "umap_cluster",
        "clustering_markers.png",
        "addmarkers-container"
      );
      document.getElementById("selectedbox").textContent =
        grabstep_umap_add_image_rslt.input_values["chosen_method"];
    }
    banned_operations();
  }
});
// banned operation
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
  document.getElementById("confirmbtn").disabled = true;
  document.getElementById("clustering_method").disabled = true;
  document.getElementById("n_neighbors").disabled = true;
  document.getElementById("resolution").disabled = true;
  document.getElementById("processbutton").disabled = true;
  document.getElementById("Addleidens").disabled = true;
  document.getElementById("pca_markers").disabled = true;
  document.getElementById("Addedmarkers").disabled = true;
}
// show selection
function show_clustering_method(preload_clustering_results) {
  const clustering_method_ul = document.getElementById("preloadul");

  clustering_method_ul.innerHTML = "";

  for (let i = 0; i < preload_clustering_results.methods.length; i++) {
    const li = document.createElement("li");
    li.innerHTML = `
                <a onclick="select('preloadul','${preload_clustering_results.methods[i]}')"
                class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2 cursor-pointer disabled:cursor-not-allowed"
            >${preload_clustering_results.methods[i]}</a>
        `;
    clustering_method_ul.appendChild(li);
  }
}
function select(s_ul, selected) {
  document.getElementById(s_ul + "_select").textContent = selected;
  document.getElementById(s_ul + "_input").value = selected;
  if (selected == "bbknn") {
    document.getElementById("n_neighbor_block").classList.add("hidden");
    document.getElementById("n_neighbors").value = 10;
    document.getElementById("resolution").value = "";
  } else {
    document.getElementById("n_neighbor_block").classList.remove("hidden");
    document.getElementById("n_neighbors").value = "";
    document.getElementById("resolution").value = "";
  }
}
function methodtextnodefix(newText) {
  methodtextnode.nodeValue = newText;
  const method = document.getElementById("method");
  method.value = newText;
}
// preocess btn click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  const form = document.getElementById("clusterform");
  const formData = new FormData(form);
  const csrftoken = getCookie("csrftoken");
  const clusterresult = await fetchAPI("/api/clustering", formData, csrftoken);

  toggleLoading(false, "processbutton");
  show_result();
  // Preload Markers
  const preloadmerkersresult = await fetchAPI(
    "/api/preloadclustermarkers",
    formData,
    csrftoken
  );
  console.log(preloadmerkersresult.data.marker_list);
  insert_marker_options(preloadmerkersresult.data.marker_list);
  document.getElementById("nextbtn").classList.remove("hidden");
}
function show_result() {
  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("imgbox").classList.remove("hidden");
  loadImage("cluster_result", "clustering_summary.png", "summary-container");
  loadImage("cluster_result", "clustering_heatmap.png", "heatmap-container");
  loadImage("cluster_result", "clustering_leidens.png", "umap-container");
  loadImage("cluster_result", "clustering_ranking.png", "ranking-container");
}
async function addleidens() {
  toggleLoading(true, "Addleidens");
  const csrftoken = getCookie("csrftoken");
  const addleidensresult = await fetchAPI(
    "/api/addumapcluster/leidens",
    0,
    csrftoken
  );
  // add to container
  loadImage(
    "umap_cluster",
    "clustering_leidens_bysample.png",
    "addleidens-container"
  );
  toggleLoading(false, "Addleidens");
}
async function addmarkers() {
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "Addedmarkers");
  const form = document.getElementById("markerform");

  const formData = new FormData(form);

  const selectedCheckboxes = document.querySelectorAll(
    ".marker-checkbox:checked"
  );
  const selectedMarkers = Array.from(selectedCheckboxes)
    .map((checkbox) => checkbox.value)
    .filter((value) => value !== "Select All");

  selectedMarkers.forEach((marker) => {
    formData.append("markers", marker);
  });
  const addmarkersresult = await fetchAPI(
    "/api/addumapcluster/markers",
    formData,
    csrftoken
  );
  loadImage("umap_cluster", "clustering_markers.png", "addmarkers-container");
  toggleLoading(false, "Addedmarkers");
}
async function prenext(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  const m_subset_results = await fetchAPI("/api/subset/merged", 0, csrftoken);
  if (m_subset_results) {
    window.location.href = event.target.href;
  }
}
