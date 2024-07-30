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
  console.log(phenotype_result);
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
  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("imgbox").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");
}
