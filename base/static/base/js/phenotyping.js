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
document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  // grab step
  const grabstep_rslt = await grabsteps(`/getSteps/phenotyping/process/`);

  if (grabstep_rslt.message != "notfound") {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    document.getElementById("btnbox").classList.remove("hidden");
    document.getElementById("nextbtn").classList.remove("hidden");
    document.getElementById("uploaded").textContent =
      grabstep_rslt.input_values["file_name"];
    banned_operations();
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
    loadImage(
      "phenotype_result",
      "phenotyping_leidens.png",
      "p-umap-container"
    );
    loadImage(
      "phenotype_result",
      "phenotyping_ranking.png",
      "p-ranking-container"
    );
    document.getElementById("imgbox").classList.remove("hidden");
    document.getElementById("phenotype_result").textContent =
      grabstep_rslt.output_values["phenotyping_result"];

    const grabstep_preloadmarkers_rslt = await grabsteps(
      `/getSteps/phenotyping/preloadmarkers/`
    );
    console.log(grabstep_preloadmarkers_rslt);
    insert_marker_options(
      grabstep_preloadmarkers_rslt.output_values["marker_list"]
    );
    const grabstep_addphenotypes_rslt = await grabsteps(
      `/getSteps/phenotyping/addumapphenotypes/`
    );
    if (grabstep_addphenotypes_rslt.message != "notfound") {
      loadImage(
        "umap_phenotyping",
        "phenotypes_bysample.png",
        "addphenotypes-container"
      );
    }

    const grabstep_addmarkers_rslt = await grabsteps(
      `/getSteps/phenotyping/addumapmarkers/`
    );
    console.log(grabstep_addmarkers_rslt);
    if (grabstep_addmarkers_rslt.message != "notfound") {
      loadImage(
        "umap_phenotyping",
        "phenotypes_markers.png",
        "addmarkers-container"
      );
      document.getElementById("selectedbox").textContent =
        grabstep_addmarkers_rslt.input_values["chosen_method"];
    }
  }
});
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
  const fileUploadLabel = document.getElementById("file_upload_l");

  fileUploadLabel.disabled = true;
  fileUploadLabel.classList.remove("bg-sky-950");
  fileUploadLabel.classList.add("bg-sky-950/50", "cursor-not-allowed");
  document.getElementById("file_upload").disabled = true;

  document.getElementById("Addedmarkers").disabled = true;
  document.getElementById("Addphenotypes").disabled = true;
  document.getElementById("phenotype_markers").disabled = true;
}
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
    "/api/preloadphenotypemarkers",
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
  toggleLoading(true, "Addphenotypes");
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
  toggleLoading(false, "Addphenotypes");
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
  toggleLoading(false, "Addedmarkers");
}
