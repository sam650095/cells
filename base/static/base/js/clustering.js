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

  // Preload Markers
  const preloadmerkersresult = await fetchAPI(
    "/api/preloadmarkers",
    formData,
    csrftoken
  );
  console.log(preloadmerkersresult);
  insert_marker_options(preloadmerkersresult.marker_list);
  document.getElementById("nextbtn").classList.remove("hidden");
}

// marker option
function insert_marker_options(markers) {
  const dropdownList = document.querySelector("#dropdownmarkers ul");
  dropdownList.innerHTML = "";

  markers.forEach((marker, index) => {
    const li = document.createElement("li");
    li.innerHTML = `
        <div class="flex p-2 rounded hover:bg-gray-100">
          <div class="flex items-center h-5">
            <input
              id="marker-${index}"
              aria-describedby="marker-text-${index}"
              type="checkbox"
              value="${marker}"
              checked
              class="marker-checkbox w-4 h-4 text-blue-600 bg-gray-100 border-gray-300 rounded focus:ring-blue-500 focus:ring-2"
            />
          </div>
          <div class="ms-2 text-sm w-full">
            <label for="marker-${index}" class="font-bold text-gray-900">
              <div>${marker}</div>
            </label>
          </div>
        </div>
      `;
    dropdownList.appendChild(li);
  });
  logSelectedCount();
  dropdownList.addEventListener("change", function (event) {
    if (event.target.type === "checkbox") {
      if (event.target.value === "Select All") {
        handleSelectAll(event.target.checked);
      } else {
        updateSelectAllState();
      }
      logSelectedCount();
    }
  });
}

function handleSelectAll(isChecked) {
  const checkboxes = document.querySelectorAll(".marker-checkbox");
  checkboxes.forEach((checkbox) => {
    checkbox.checked = isChecked;
  });
}

function updateSelectAllState() {
  const checkboxes = document.querySelectorAll(
    '.marker-checkbox:not([value="Select All"])'
  );
  const selectAllCheckbox = document.querySelector(
    '.marker-checkbox[value="Select All"]'
  );
  const allChecked = Array.from(checkboxes).every(
    (checkbox) => checkbox.checked
  );
  selectAllCheckbox.checked = allChecked;
}

function logSelectedCount() {
  const selectedCheckboxes = document.querySelectorAll(
    ".marker-checkbox:checked"
  );
  const selectedMarkers = Array.from(selectedCheckboxes)
    .map((checkbox) => checkbox.value)
    .filter((value) => value !== "Select All");

  const selectedCount = selectedMarkers.length;

  const selectedbox = document.getElementById("selectedbox");
  selectedbox.textContent = `Currently selected: ${selectedCount} markers`;
  const markersdiv = document.createElement("div");
  markersdiv.textContent = `[${selectedMarkers.join(", ")}]`;
  markersdiv.classList.add("indent-8");

  selectedbox.appendChild(markersdiv);
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
