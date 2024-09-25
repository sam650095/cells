document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const markerlist_results = await fetchAPI("/api/preloadpca", 0, csrftoken);
  insert_marker_options(markerlist_results.data.marker_list);
  // grab step
  const grabstep_rslt = await grabsteps(`/getSteps/pca/process/`);

  if (grabstep_rslt.message != "notfound") {
    stepped = true;
    document.getElementById("watchonly").classList.remove("hidden");
    grabstep_rslt.input_values["chosen_markers"].forEach((element) => {
      document.getElementById("marker-" + element).checked = true;
    });
    logSelectedCount();
    updateSelectAllState();
    show_list(grabstep_rslt.output_values);
    banned_operations();
  }
});
function banned_operations() {
  const processBtn = document.getElementById("processbutton");
  processBtn.disabled = true;
  const pca_markers = document.getElementById("pca_markers");
  pca_markers.disabled = true;
}
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
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

  const pca_box = document.getElementById("pca_box");
  pca_box.textContent = "";
  const process_result = await fetchAPI("/api/pca", formData, csrftoken);
  console.log(process_result);
  show_list(process_result.data);
}
function show_list(process_result) {
  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-fit"
  );
  const mli = document.createElement("li");
  mli.textContent = process_result["merged_results"];
  const nli = document.createElement("li");
  nli.textContent = process_result["n_pcs_results"];
  const img = document.createElement("img");
  img.src = `/media/pca_img/${process_result["save_img_name"]}?t=${Date.now()}`;
  img.classList.add("m-5");
  nli.appendChild(img);
  ul.appendChild(mli);
  ul.appendChild(nli);
  pca_box.appendChild(ul);

  // next page btn
  toggleLoading(false, "processbutton");
  document.getElementById("nextbtn").classList.remove("hidden");
}
