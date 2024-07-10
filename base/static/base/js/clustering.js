let methodtextnode;
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
  toggleLoading(true);
  const form = document.getElementById("clusterform");
  const formData = new FormData(form);
  console.log(formData);
  const csrftoken = getCookie("csrftoken");
  const clusterresult = await fetchAPI("/api/clustering", formData, csrftoken);
  console.log(clusterresult);

  toggleLoading(false);
  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("imgbox").classList.remove("hidden");
  loadImage("cluster_result", "clustering_summary.png", "summary-container");
  loadImage("cluster_result", "clustering_heatmap.png", "heatmap-container");
  loadImage("cluster_result", "clustering_leidens.png", "umap-container");
  loadImage("cluster_result", "clustering_ranking.png", "ranking-container");
}
