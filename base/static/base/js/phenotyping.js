const fileInput = document.getElementById("file_upload");
const uploadedDiv = document.getElementById("uploaded");
let selectedFile = null;
fileInput.addEventListener("change", function (event) {
  const file = event.target.files[0];
  if (file) {
    uploadedDiv.textContent = file.name;
    selectedFile = file;
  } else {
    uploadedDiv.textContent = "";
    selectedFile = null;
  }
});
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const formData = new FormData();
  formData.append("file", selectedFile);

  const phenotype_result = await fetchAPI(
    "/api/phenotyping",
    formData,
    csrftoken
  );
  document.getElementById("phenotype_result").textContent =
    phenotype_result.data.phenotyping_result;
  toggleLoading(false, "processbutton");
  loadImage(
    "phenotype_result",
    "phenotyping_summary.png",
    "p-summary-container"
  );
  loadImage(
    "phenotype_result",
    "phenotyping_heatmap.png",
    "p-heatmap-container"
  );
  loadImage("phenotype_result", "phenotyping_leidens.png", "p-umap-container");
  loadImage(
    "phenotype_result",
    "phenotyping_ranking.png",
    "p-ranking-container"
  );
  // Preload Markers
  const preloadmerkersresult = await fetchAPI(
    "/api/preloadclustermarkers",
    formData,
    csrftoken
  );
  console.log(preloadmerkersresult);
  // Insert markers options
  insert_marker_options(preloadmerkersresult.data.marker_list);
  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("imgbox").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");
}

async function addphenotypes() {
  const csrftoken = getCookie("csrftoken");
  const addphenotypesresult = await fetchAPI(
    "/api/addumapphenotype/phenotypes",
    0,
    csrftoken
  );
  console.log(addphenotypesresult);
  // add to container
  loadImage(
    "umap_phenotyping",
    "phenotypes_bysample.png",
    "addphenotypes-container"
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
    "/api/addumapphenotype/markers",
    formData,
    csrftoken
  );
  console.log(addleidensresult);
  loadImage(
    "umap_phenotyping",
    "phenotypes_markers.png",
    "addmarkers-container"
  );
}
