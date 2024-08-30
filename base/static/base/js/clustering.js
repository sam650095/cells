let methodtextnode;
async function prenext(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  const m_subset_results = await fetchAPI("/api/subset/merged", 0, csrftoken);
  console.log(m_subset_results);
  if (m_subset_results) {
    window.location.href = event.target.href;
  }
}
document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_clustering_results = await fetchAPI(
    "/api/preloadclustering",
    0,
    csrftoken
  );
  console.log(preload_clustering_results);
  show_clustering_method(preload_clustering_results);
});
// show selection
function show_clustering_method(preload_clustering_results) {
  const clustering_method = document.getElementById("clustering_method");
  const clustering_method_ul = document.getElementById("preloadul");

  clustering_method_ul.innerHTML = "";
  methodtextnode = document.createTextNode("none");
  clustering_method.insertBefore(methodtextnode, clustering_method.firstChild);
  for (let i = 0; i < preload_clustering_results.methods.length; i++) {
    const li = document.createElement("li");
    li.innerHTML = `
                <a onclick="methodtextnodefix('${preload_clustering_results.methods[i]}')"
                class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
            >${preload_clustering_results.methods[i]}</a>
        `;
    clustering_method_ul.appendChild(li);
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
  console.log(formData);
  const csrftoken = getCookie("csrftoken");
  const clusterresult = await fetchAPI("/api/clustering", formData, csrftoken);
  console.log(clusterresult);

  toggleLoading(false, "processbutton");
  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("imgbox").classList.remove("hidden");
  loadImage("cluster_result", "clustering_summary.png", "summary-container");
  loadImage("cluster_result", "clustering_heatmap.png", "heatmap-container");
  loadImage("cluster_result", "clustering_leidens.png", "umap-container");
  loadImage("cluster_result", "clustering_ranking.png", "ranking-container");

  // Preload Markers
  const preloadmerkersresult = await fetchAPI(
    "/api/preloadclustermarkers",
    formData,
    csrftoken
  );
  console.log(preloadmerkersresult);
  insert_marker_options(preloadmerkersresult.marker_list);
  document.getElementById("nextbtn").classList.remove("hidden");
}

async function addleidens() {
  const csrftoken = getCookie("csrftoken");
  const addleidensresult = await fetchAPI(
    "/api/addumapcluster/leidens",
    0,
    csrftoken
  );
  console.log(addleidensresult);
  // add to container
  loadImage(
    "umap_cluster",
    "clustering_leidens_bysample.png",
    "addleidens-container"
  );
}
async function addmarkers() {
  const csrftoken = getCookie("csrftoken");

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
  const addleidensresult = await fetchAPI(
    "/api/addumapcluster/markers",
    formData,
    csrftoken
  );
  console.log(addleidensresult);
  loadImage("umap_cluster", "clustering_markers.png", "addmarkers-container");
}
