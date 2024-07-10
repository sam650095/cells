document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const markerlist_results = await fetchAPI("/api/preloadpca", 0, csrftoken);
  console.log(markerlist_results);
  insert_marker_options(markerlist_results.marker_list);
});

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
  console.log(`Currently selected: ${selectedCount} marker(s)`);
  console.log("Selected markers:", selectedMarkers);

  const selectedbox = document.getElementById("selectedbox");
  selectedbox.textContent = "";
  const ul = document.createElement("ul");
  ul.classList.add(
    "space-y-2",
    "text-gray-500",
    "list-disc",
    "list-inside",
    "w-full",
    "break-words"
  );

  const li = document.createElement("li");
  li.textContent = `Currently selected: ${selectedCount} markers`;
  const markersdiv = document.createElement("div");
  markersdiv.textContent = `[${selectedMarkers.join(", ")}]`;
  markersdiv.classList.add("indent-8");
  li.appendChild(markersdiv);
  ul.appendChild(li);

  selectedbox.appendChild(ul);
}
// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true);
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
  toggleLoading(false);
  document.getElementById("nextbtn").classList.remove("hidden");
}
